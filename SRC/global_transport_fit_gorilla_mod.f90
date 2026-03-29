module global_transport_fit_gorilla_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    public :: run_gorilla_mc_global_transport_fit

contains

subroutine run_gorilla_mc_global_transport_fit(control, result)

    use field_mod, only: ipert
    use global_transport_fit_core_mod, only: fit_global_transport
    use global_transport_fit_io_mod, only: write_global_transport_fit_outputs
    use global_transport_fit_math_mod, only: compute_geometry_from_boundaries
    use global_transport_fit_settings_mod, only: filename_comparison_profiles, filename_experiment_1_summary, &
        filename_experiment_2_summary, filename_fit_profiles, filename_fit_summary, filename_local_profiles, &
        local_density_gradient_1, local_density_gradient_2, local_density_profile_reference_s, local_v_E, &
        n_local_reference_surfaces, source_gradient_1, source_gradient_2, source_reference_s_1, source_reference_s_2
    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_fit_control_t, &
        global_transport_fit_result_t
    use gorilla_applets_settings_mod, only: i_option
    use gorilla_applets_types_mod, only: in, s
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use transport_benchmark_utils_mod, only: prepare_minimal_transport_output
    use utils_data_pre_and_post_processing_mod, only: get_ipert, set_seed_for_random_numbers
    use utils_self_consistent_ef_mod, only: allocate_electric_potential_type, associate_flux_labels_with_tetrahedra_and_vertices, &
        calc_s_shell_volumes, print_errors_for_bad_inputs, read_self_consistent_electric_field_inp_into_type
    use volume_integrals_and_sqrt_g_mod, only: calc_volume_integrals_in_flux_coordinates

    type(global_transport_fit_control_t), intent(in) :: control
    type(global_transport_experiment_t), allocatable :: experiments(:)
    type(global_transport_fit_result_t), intent(out) :: result
    real(dp), allocatable :: boundary_areas(:)
    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: cell_centers(:)
    real(dp), allocatable :: local_a(:)
    real(dp), allocatable :: local_b(:)
    real(dp), allocatable :: local_boundary_s(:)
    real(dp), allocatable :: shell_volumes(:)
    integer, allocatable :: local_boundary_indices(:)

    call set_seed_for_random_numbers()
    call read_self_consistent_electric_field_inp_into_type()
    call get_ipert()
    call initialize_gorilla(i_option, ipert)

    in%tracer_species = 1
    in%n_species = 1
    in%boole_custom_background = .true.
    in%background_density_cm3 = (/in%density, in%density/)
    in%background_temperature_eV = (/in%energy_eV, in%energy_eV/)
    in%collision_operator = 4
    s%temperature = in%energy_eV
    call prepare_minimal_transport_output()
    call calc_volume_integrals_in_flux_coordinates()
    call associate_flux_labels_with_tetrahedra_and_vertices()
    call allocate_electric_potential_type()
    call calc_s_shell_volumes()
    call print_errors_for_bad_inputs()

    call extract_boundary_geometry(boundary_s, shell_volumes)
    call compute_geometry_from_boundaries(boundary_s, shell_volumes, cell_centers, boundary_areas)

    allocate(experiments(2))
    call run_electron_source_experiment(source_gradient_1, source_reference_s_1, boundary_s, shell_volumes, boundary_areas, &
        filename_experiment_1_summary, experiments(1))
    call run_electron_source_experiment(source_gradient_2, source_reference_s_2, boundary_s, shell_volumes, boundary_areas, &
        filename_experiment_2_summary, experiments(2))

    call fit_global_transport(experiments, control, result)
    call write_global_transport_fit_outputs(boundary_s, result, filename_fit_summary, filename_fit_profiles)

    call build_local_transport_reference(boundary_s, n_local_reference_surfaces, local_density_gradient_1, &
        local_density_gradient_2, local_density_profile_reference_s, local_v_E, local_boundary_indices, local_boundary_s, &
        local_a, local_b)
    call write_local_transport_profiles(local_boundary_indices, local_boundary_s, local_a, local_b, filename_local_profiles)
    call write_transport_comparison(local_boundary_indices, local_boundary_s, local_a, local_b, result%a_profile, result%a_std, &
        result%b_profile, result%b_std, filename_comparison_profiles)

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
    summary_filename, experiment)

    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_samples_t
    use gorilla_applets_types_mod, only: in, output
    use tetra_grid_settings_mod, only: grid_size
    use utils_data_pre_and_post_processing_mod, only: calc_collision_coefficients_for_all_tetrahedra, initialize_exit_data, &
        prepare_next_round_of_parallelised_particle_pushing
    use utils_self_consistent_ef_mod, only: allocate_start_type, parallelised_particle_pushing, &
        set_particle_type_specifications, set_rest_of_individual_particle_specifications, set_starting_positions

    real(dp), intent(in) :: source_gradient
    real(dp), intent(in) :: source_reference_s
    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: shell_volumes(:)
    character(len=*), intent(in) :: summary_filename
    type(global_transport_experiment_t), intent(out) :: experiment

    type(global_transport_samples_t) :: samples
    real(dp), allocatable :: rand_matrix(:, :, :)

    allocate(samples%source_weight_sum(grid_size(1)))
    allocate(samples%shell_time_sum(grid_size(1)))
    allocate(samples%shell_time_sumsq(grid_size(1)))
    allocate(samples%boundary_weighted_flux_sum(grid_size(1) + 1))
    allocate(samples%boundary_weighted_flux_sumsq(grid_size(1) + 1))
    samples%source_weight_sum = 0.0_dp
    samples%shell_time_sum = 0.0_dp
    samples%shell_time_sumsq = 0.0_dp
    samples%boundary_weighted_flux_sum = 0.0_dp
    samples%boundary_weighted_flux_sumsq = 0.0_dp

    allocate(rand_matrix(5, in%num_particles, 1))
    call random_number(rand_matrix)

    output%radial_flux = 0.0_dp
    output%weighted_radial_flux = 0.0_dp

    in%density_log_gradient_per_s = source_gradient
    in%density_profile_reference_s = source_reference_s

    call allocate_start_type(in%num_particles)
    call set_particle_type_specifications()
    call initialize_exit_data(in%num_particles)
    call set_starting_positions(rand_matrix, species_in=(/1/))
    call set_rest_of_individual_particle_specifications(rand_matrix, boole_diffusion_coefficient_in=.true., &
        species_in=(/1/), n_particles_in=in%num_particles)
    call prepare_next_round_of_parallelised_particle_pushing(1)
    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra(1)
    call parallelised_particle_pushing(species=1, j=1, boole_diffusion_coefficient=.false., n_particles_in=in%num_particles, &
        transport_samples=samples)

    call build_experiment_from_samples(boundary_s, shell_volumes, boundary_areas, samples, output%weighted_radial_flux, experiment)
    call write_experiment_summary(boundary_s, shell_volumes, boundary_areas, samples, output%weighted_radial_flux, &
        summary_filename, source_gradient, source_reference_s)

end subroutine run_electron_source_experiment

subroutine build_experiment_from_samples(boundary_s, shell_volumes, boundary_areas, samples, boundary_flux_sum, experiment)

    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_samples_t

    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: boundary_flux_sum(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: shell_volumes(:)
    type(global_transport_samples_t), intent(in) :: samples
    type(global_transport_experiment_t), intent(out) :: experiment

    real(dp) :: normalization
    real(dp) :: particles

    particles = real(max(samples%n_particles, 1), dp)
    normalization = max(sum(samples%source_weight_sum), 1.0d-20)

    experiment%boundary_s = boundary_s
    experiment%shell_volumes = shell_volumes
    experiment%source = samples%source_weight_sum / (normalization * shell_volumes)
    experiment%density = samples%shell_time_sum / (normalization * shell_volumes)
    allocate(experiment%density_variance(size(shell_volumes)))
    experiment%density_variance = max(samples%shell_time_sumsq / particles - (samples%shell_time_sum / particles)**2, 0.0_dp)
    experiment%density_variance = experiment%density_variance / &
        (particles * (normalization / particles)**2 * shell_volumes**2)
    experiment%flux = boundary_flux_sum / (normalization * boundary_areas)
    allocate(experiment%flux_variance(size(boundary_s)))
    experiment%flux_variance = max(samples%boundary_weighted_flux_sumsq / particles - &
        (samples%boundary_weighted_flux_sum / particles)**2, 0.0_dp)
    experiment%flux_variance = experiment%flux_variance / &
        (particles * (normalization / particles)**2 * boundary_areas**2)
    experiment%flux_variance = max(experiment%flux_variance, 1.0d-30)

end subroutine build_experiment_from_samples

subroutine write_experiment_summary(boundary_s, shell_volumes, boundary_areas, samples, boundary_flux_sum, filename, &
    source_gradient, source_reference_s)

    use global_transport_fit_types_mod, only: global_transport_samples_t

    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: boundary_flux_sum(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: shell_volumes(:)
    character(len=*), intent(in) :: filename
    type(global_transport_samples_t), intent(in) :: samples
    real(dp), intent(in) :: source_gradient
    real(dp), intent(in) :: source_reference_s

    integer :: i
    integer :: io_unit
    real(dp) :: normalization

    normalization = max(sum(samples%source_weight_sum), 1.0d-20)
    open(newunit=io_unit, file=filename, status='replace', action='write')
    write(io_unit, '(A)') 'shell_index shell_center_s source_density density boundary_flux_density'
    do i = 1, size(shell_volumes)
        write(io_unit, '(I0,4(1X,ES24.16))') i, 0.5_dp * (boundary_s(i) + boundary_s(i + 1)), &
            samples%source_weight_sum(i) / (normalization * shell_volumes(i)), &
            samples%shell_time_sum(i) / (normalization * shell_volumes(i)), &
            boundary_flux_sum(i + 1) / (normalization * boundary_areas(i + 1))
    end do
    write(io_unit, '(A)') 'source_gradient source_reference_s total_source_weight lost_particles n_particles'
    write(io_unit, '(5(1X,ES24.16))') source_gradient, source_reference_s, normalization, real(samples%lost_particles, dp), &
        real(samples%n_particles, dp)
    close(io_unit)

end subroutine write_experiment_summary

subroutine build_local_transport_reference(boundary_s, n_reference_surfaces, density_gradient_1, density_gradient_2, &
    density_profile_reference_s, v_E, boundary_indices, local_boundary_s, local_a, local_b)

    use gorilla_applets_types_mod, only: in
    use transport_benchmark_utils_mod, only: recover_transport_coefficients_from_flux_pair
    use usual_transport_benchmark_mod, only: run_usual_transport_benchmark_case, usual_transport_result_t

    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: density_gradient_1
    real(dp), intent(in) :: density_gradient_2
    real(dp), intent(in) :: density_profile_reference_s
    real(dp), intent(in) :: v_E
    integer, intent(in) :: n_reference_surfaces
    integer, allocatable, intent(out) :: boundary_indices(:)
    real(dp), allocatable, intent(out) :: local_boundary_s(:)
    real(dp), allocatable, intent(out) :: local_a(:)
    real(dp), allocatable, intent(out) :: local_b(:)

    type(usual_transport_result_t) :: result_1
    type(usual_transport_result_t) :: result_2
    real(dp) :: density_1
    real(dp) :: density_2
    integer :: i

    call select_local_reference_boundaries(boundary_s, n_reference_surfaces, boundary_indices)
    allocate(local_boundary_s(size(boundary_indices)))
    allocate(local_a(size(boundary_indices)))
    allocate(local_b(size(boundary_indices)))

    do i = 1, size(boundary_indices)
        call run_usual_transport_benchmark_case(boundary_s(boundary_indices(i)), 1, in%density, in%density, in%energy_eV, &
            in%density, in%energy_eV, in%energy_eV, in%time_step, in%num_particles, density_gradient_1, &
            density_profile_reference_s, v_E, in%collision_operator, in%boole_precalc_collisions, in%i_integrator_type, &
            result_1)
        call run_usual_transport_benchmark_case(boundary_s(boundary_indices(i)), 1, in%density, in%density, in%energy_eV, &
            in%density, in%energy_eV, in%energy_eV, in%time_step, in%num_particles, density_gradient_2, &
            density_profile_reference_s, v_E, in%collision_operator, in%boole_precalc_collisions, in%i_integrator_type, &
            result_2)
        if (abs(result_1%target_boundary_s - result_2%target_boundary_s) > 1.0d-12) then
            print *, 'Error: local transport reference runs selected different target boundaries.'
            stop
        end if
        density_1 = in%density * exp(density_gradient_1 * (result_1%target_boundary_s - density_profile_reference_s))
        density_2 = in%density * exp(density_gradient_2 * (result_2%target_boundary_s - density_profile_reference_s))
        call recover_transport_coefficients_from_flux_pair(result_1%weighted_flux_density, result_2%weighted_flux_density, &
            density_1, density_2, density_gradient_1, density_gradient_2, local_a(i), local_b(i))
        local_boundary_s(i) = result_1%target_boundary_s
    end do

end subroutine build_local_transport_reference

subroutine select_local_reference_boundaries(boundary_s, n_reference_surfaces, boundary_indices)

    real(dp), intent(in) :: boundary_s(:)
    integer, intent(in) :: n_reference_surfaces
    integer, allocatable, intent(out) :: boundary_indices(:)

    integer :: i
    integer :: n_interior
    integer :: n_selected

    n_interior = max(size(boundary_s) - 2, 1)
    n_selected = min(max(n_reference_surfaces, 1), n_interior)
    allocate(boundary_indices(n_selected))

    if (n_selected == 1) then
        boundary_indices(1) = 1 + (n_interior + 1) / 2
        return
    end if

    do i = 1, n_selected
        boundary_indices(i) = 2 + int(real(i - 1, dp) * real(n_interior - 1, dp) / real(n_selected - 1, dp))
    end do

end subroutine select_local_reference_boundaries

subroutine write_local_transport_profiles(boundary_indices, boundary_s, local_a, local_b, filename)

    integer, intent(in) :: boundary_indices(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: local_a(:)
    real(dp), intent(in) :: local_b(:)
    character(len=*), intent(in) :: filename

    integer :: i
    integer :: io_unit

    open(newunit=io_unit, file=filename, status='replace', action='write')
    write(io_unit, '(A)') 'boundary_index boundary_s local_A local_B'
    do i = 1, size(boundary_indices)
        write(io_unit, '(I0,3(1X,ES24.16))') boundary_indices(i), boundary_s(i), local_a(i), local_b(i)
    end do
    close(io_unit)

end subroutine write_local_transport_profiles

subroutine write_transport_comparison(boundary_indices, boundary_s, local_a, local_b, fit_a, fit_a_std, fit_b, fit_b_std, &
    filename)

    integer, intent(in) :: boundary_indices(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: local_a(:)
    real(dp), intent(in) :: local_b(:)
    real(dp), intent(in) :: fit_a(:)
    real(dp), intent(in) :: fit_a_std(:)
    real(dp), intent(in) :: fit_b(:)
    real(dp), intent(in) :: fit_b_std(:)
    character(len=*), intent(in) :: filename

    integer :: i
    integer :: io_unit

    open(newunit=io_unit, file=filename, status='replace', action='write')
    write(io_unit, '(A)') 'boundary_index boundary_s local_A fit_A fit_A_std local_B fit_B fit_B_std'
    do i = 1, size(boundary_indices)
        write(io_unit, '(I0,7(1X,ES24.16))') boundary_indices(i), boundary_s(i), local_a(i), fit_a(boundary_indices(i)), &
            fit_a_std(boundary_indices(i)), local_b(i), fit_b(boundary_indices(i)), fit_b_std(boundary_indices(i))
    end do
    close(io_unit)

end subroutine write_transport_comparison

end module global_transport_fit_gorilla_mod
