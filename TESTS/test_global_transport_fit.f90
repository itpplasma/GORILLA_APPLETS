program test_global_transport_fit

    use, intrinsic :: iso_fortran_env, only: dp => real64

    use global_transport_fit_core_mod, only: compute_objective_gradient_and_jacobian, fit_global_transport, &
        generate_synthetic_experiment
    use global_transport_fit_math_mod, only: build_piecewise_linear_basis, compute_geometry_from_boundaries
    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_fit_control_t, &
        global_transport_fit_result_t
    use transport_benchmark_utils_mod, only: recover_transport_coefficients_from_flux_pair
    use transport_statistics_mod, only: build_supported_boundary_mask, compute_mean_and_variance, &
        fit_linear_transport_response, sanitize_transport_output_triplet, transport_signal_supported

    implicit none

    call test_fit_recovers_manufactured_profiles()
    call test_fit_recovers_manufactured_profiles_density_only()
    call test_gradient_matches_finite_difference()
    call test_forward_problem_has_zero_inner_flux()
    call test_boundary_geometry_has_zero_axis_area()
    call test_supported_boundary_mask_keeps_outer_loss_boundary()
    call test_transport_signal_requires_density_separation()
    call test_transport_output_sanitization_handles_invalid_sigma()
    call test_local_flux_pair_recovery()
    call test_local_reference_boundaries_avoid_axis_and_edge()
    call test_weighted_linear_transport_response()
    call test_batch_mean_and_variance()
    call test_convergence_history_populated()

contains

subroutine test_fit_recovers_manufactured_profiles()

    type(global_transport_experiment_t), allocatable :: experiments(:)
    type(global_transport_fit_control_t) :: control
    type(global_transport_fit_result_t) :: result
    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: shell_volumes(:)
    real(dp), allocatable :: source_1(:)
    real(dp), allocatable :: source_2(:)
    real(dp), allocatable :: knot_s(:)
    real(dp), allocatable :: basis(:, :)
    real(dp), allocatable :: true_a(:)
    real(dp), allocatable :: true_b(:)
    real(dp), allocatable :: expected_a(:)
    real(dp), allocatable :: expected_b(:)

    call make_test_geometry(boundary_s, shell_volumes, source_1, source_2)

    control%n_knots_a = 3
    control%n_knots_b = 3
    control%max_lm_iterations = 25
    control%regularization_a = 1.0d-12
    control%regularization_b = 1.0d-12
    control%lm_damping = 1.0d-6

    call make_knot_points(boundary_s, control%n_knots_a, knot_s)
    allocate(true_a(control%n_knots_a))
    allocate(true_b(control%n_knots_b))
    true_a = (/5.0d-2, 2.0d-2, -1.0d-2/)
    true_b = log((/2.0d-2, 3.0d-2, 5.0d-2/))

    allocate(experiments(2))
    call generate_synthetic_experiment(boundary_s, shell_volumes, source_1, knot_s, knot_s, true_a, true_b, 0.0_dp, experiments(1))
    call generate_synthetic_experiment(boundary_s, shell_volumes, source_2, knot_s, knot_s, true_a, true_b, 0.0_dp, experiments(2))

    call fit_global_transport(experiments, control, result)

    call build_piecewise_linear_basis(boundary_s, knot_s, basis)
    expected_a = matmul(basis, true_a)
    expected_b = exp(matmul(basis, true_b))

    if (.not. result%converged) then
        print *, 'objective = ', result%objective
        print *, 'iterations = ', result%n_iterations
        print *, 'gradient = ', result%gradient
        print *, 'A profile = ', result%a_profile
        print *, 'B profile = ', result%b_profile
    end if
    call assert_true(result%converged, 'LM fit did not converge')
    call assert_close(maxval(abs(result%a_profile - expected_a)), 0.0_dp, 1.0d-6, 'A profile recovery failed')
    call assert_close(maxval(abs(result%b_profile - expected_b)), 0.0_dp, 1.0d-6, 'B profile recovery failed')
    call assert_true(all(result%a_std >= 0.0_dp), 'A profile standard deviations must be non-negative')
    call assert_true(all(result%b_std >= 0.0_dp), 'B profile standard deviations must be non-negative')

end subroutine test_fit_recovers_manufactured_profiles

subroutine test_fit_recovers_manufactured_profiles_density_only()

    type(global_transport_experiment_t), allocatable :: experiments(:)
    type(global_transport_fit_control_t) :: control
    type(global_transport_fit_result_t) :: result
    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: shell_volumes(:)
    real(dp), allocatable :: source_1(:)
    real(dp), allocatable :: source_2(:)
    real(dp), allocatable :: knot_s(:)
    real(dp), allocatable :: basis(:, :)
    real(dp), allocatable :: true_a(:)
    real(dp), allocatable :: true_b(:)
    real(dp), allocatable :: expected_a(:)
    real(dp), allocatable :: expected_b(:)

    call make_test_geometry(boundary_s, shell_volumes, source_1, source_2)

    control%n_knots_a = 3
    control%n_knots_b = 3
    control%max_lm_iterations = 25
    control%regularization_a = 1.0d-12
    control%regularization_b = 1.0d-12
    control%lm_damping = 1.0d-6
    control%use_flux_objective = .false.

    call make_knot_points(boundary_s, control%n_knots_a, knot_s)
    allocate(true_a(control%n_knots_a))
    allocate(true_b(control%n_knots_b))
    true_a = (/5.0d-2, 2.0d-2, -1.0d-2/)
    true_b = log((/2.0d-2, 3.0d-2, 5.0d-2/))

    allocate(experiments(2))
    call generate_synthetic_experiment(boundary_s, shell_volumes, source_1, knot_s, knot_s, true_a, true_b, 0.0_dp, experiments(1))
    call generate_synthetic_experiment(boundary_s, shell_volumes, source_2, knot_s, knot_s, true_a, true_b, 0.0_dp, experiments(2))

    call fit_global_transport(experiments, control, result)

    call build_piecewise_linear_basis(boundary_s, knot_s, basis)
    expected_a = matmul(basis, true_a)
    expected_b = exp(matmul(basis, true_b))

    call assert_true(result%converged, 'Density-only LM fit did not converge')
    call assert_close(maxval(abs(result%a_profile - expected_a)), 0.0_dp, 1.0d-6, 'Density-only A profile recovery failed')
    call assert_close(maxval(abs(result%b_profile - expected_b)), 0.0_dp, 1.0d-6, 'Density-only B profile recovery failed')

end subroutine test_fit_recovers_manufactured_profiles_density_only

subroutine test_gradient_matches_finite_difference()

    type(global_transport_experiment_t), allocatable :: experiments(:)
    type(global_transport_fit_control_t) :: control
    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: shell_volumes(:)
    real(dp), allocatable :: source_1(:)
    real(dp), allocatable :: source_2(:)
    real(dp), allocatable :: knot_s(:)
    real(dp), allocatable :: basis(:, :)
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: true_a(:)
    real(dp), allocatable :: true_b(:)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: jacobian(:, :)
    real(dp) :: objective
    real(dp) :: epsilon
    real(dp) :: fd_objective_minus
    real(dp) :: fd_objective_plus
    real(dp), allocatable :: fd_gradient(:)
    integer :: i

    call make_test_geometry(boundary_s, shell_volumes, source_1, source_2)

    control%n_knots_a = 3
    control%n_knots_b = 3
    control%regularization_a = 1.0d-10
    control%regularization_b = 1.0d-10

    call make_knot_points(boundary_s, control%n_knots_a, knot_s)
    call build_piecewise_linear_basis(boundary_s, knot_s, basis)
    allocate(true_a(3))
    allocate(true_b(3))
    true_a = (/5.0d-2, 2.0d-2, -1.0d-2/)
    true_b = log((/2.0d-2, 3.0d-2, 5.0d-2/))

    allocate(experiments(2))
    call generate_synthetic_experiment(boundary_s, shell_volumes, source_1, knot_s, knot_s, true_a, true_b, 0.0_dp, experiments(1))
    call generate_synthetic_experiment(boundary_s, shell_volumes, source_2, knot_s, knot_s, true_a, true_b, 0.0_dp, experiments(2))

    allocate(parameters(6))
    parameters = (/2.0d-2, -1.0d-2, 1.0d-2, log(1.5d-2), log(2.5d-2), log(4.0d-2)/)
    call compute_objective_gradient_and_jacobian(experiments, control, basis, basis, parameters, objective, gradient, jacobian)

    epsilon = 1.0d-7
    allocate(fd_gradient(size(parameters)))
    do i = 1, size(parameters)
        parameters(i) = parameters(i) + epsilon
        call compute_objective_only(experiments, control, basis, parameters, fd_objective_plus)
        parameters(i) = parameters(i) - 2.0_dp * epsilon
        call compute_objective_only(experiments, control, basis, parameters, fd_objective_minus)
        parameters(i) = parameters(i) + epsilon
        fd_gradient(i) = (fd_objective_plus - fd_objective_minus) / (2.0_dp * epsilon)
    end do

    call assert_relative_close(maxval(abs(gradient - fd_gradient)), maxval(abs(fd_gradient)), 1.0d-9, 'LM gradient mismatch')

end subroutine test_gradient_matches_finite_difference

subroutine test_forward_problem_has_zero_inner_flux()

    type(global_transport_experiment_t) :: experiment
    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: shell_volumes(:)
    real(dp), allocatable :: source_1(:)
    real(dp), allocatable :: source_2(:)
    real(dp), allocatable :: knot_s(:)
    real(dp), allocatable :: true_a(:)
    real(dp), allocatable :: true_b(:)

    call make_test_geometry(boundary_s, shell_volumes, source_1, source_2)
    call make_knot_points(boundary_s, 3, knot_s)
    allocate(true_a(3))
    allocate(true_b(3))
    true_a = (/5.0d-2, 2.0d-2, -1.0d-2/)
    true_b = log((/2.0d-2, 3.0d-2, 5.0d-2/))

    call generate_synthetic_experiment(boundary_s, shell_volumes, source_1, knot_s, knot_s, true_a, true_b, 0.0_dp, experiment)

    call assert_close(experiment%flux(1), 0.0_dp, 1.0d-12, 'Inner boundary flux must vanish for torus topology')

end subroutine test_forward_problem_has_zero_inner_flux

subroutine test_boundary_geometry_has_zero_axis_area()

    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: boundary_areas(:)
    real(dp), allocatable :: cell_centers(:)
    real(dp), allocatable :: cell_volume_derivative(:)
    real(dp), allocatable :: shell_volumes(:)
    real(dp), allocatable :: source_1(:)
    real(dp), allocatable :: source_2(:)

    integer :: i

    call make_test_geometry(boundary_s, shell_volumes, source_1, source_2)
    call compute_geometry_from_boundaries(boundary_s, shell_volumes, cell_centers, boundary_areas)

    call assert_close(boundary_areas(1), 0.0_dp, 1.0d-12, 'Magnetic-axis boundary area must vanish')
    call assert_true(all(boundary_areas(2:) > 0.0_dp), 'All non-axis transport geometry factors must stay positive')

    allocate(cell_volume_derivative(size(shell_volumes)))
    do i = 1, size(shell_volumes)
        cell_volume_derivative(i) = shell_volumes(i) / (boundary_s(i + 1) - boundary_s(i))
    end do

    call assert_close(boundary_areas(size(boundary_areas)), cell_volume_derivative(size(cell_volume_derivative)), 1.0d-12, &
        'Outer transport geometry factor must match the outer shell volume derivative')
    do i = 2, size(boundary_areas) - 1
        call assert_close(boundary_areas(i), 0.5_dp * (cell_volume_derivative(i - 1) + cell_volume_derivative(i)), 1.0d-12, &
            'Interior transport geometry factor must average adjacent shell volume derivatives')
    end do

end subroutine test_boundary_geometry_has_zero_axis_area

subroutine test_supported_boundary_mask_keeps_outer_loss_boundary()

    real(dp), allocatable :: flux(:)
    real(dp), allocatable :: flux_variance(:)
    logical, allocatable :: supported_mask(:)

    flux = (/0.0_dp, 5.0d-4, 0.0_dp, 2.0d-3/)
    flux_variance = (/0.0_dp, 1.0d-8, 1.0d-8, 1.0d-8/)

    call build_supported_boundary_mask(flux, flux_variance, 2.0_dp, supported_mask)

    call assert_true(.not. supported_mask(1), 'Magnetic-axis boundary must stay unsupported')
    call assert_true(supported_mask(2), 'Interior supported boundary was dropped')
    call assert_true(.not. supported_mask(3), 'Noise-only boundary must stay unsupported')
    call assert_true(supported_mask(4), 'Outer absorbing-loss boundary must stay supported')

end subroutine test_supported_boundary_mask_keeps_outer_loss_boundary

subroutine test_transport_signal_requires_density_separation()

    call assert_true(.not. transport_signal_supported(3, 2, 4.0d-2, 5.0d-2), &
        'Source-dominated experiment must not pass the transport-signal gate')
    call assert_true(transport_signal_supported(3, 2, 6.0d-2, 5.0d-2), &
        'Transport-informative experiment must pass the transport-signal gate')

end subroutine test_transport_signal_requires_density_separation

subroutine test_transport_output_sanitization_handles_invalid_sigma()

    real(dp) :: sanitized_2sigma
    real(dp) :: sanitized_mean
    real(dp) :: sanitized_std

    call sanitize_transport_output_triplet(.false., 1.0_dp, huge(1.0_dp), sanitized_mean, sanitized_std, sanitized_2sigma)

    call assert_true(sanitized_mean /= sanitized_mean, 'Invalid transport mean must be written as NaN')
    call assert_true(sanitized_std /= sanitized_std, 'Invalid transport standard deviation must be written as NaN')
    call assert_true(sanitized_2sigma /= sanitized_2sigma, 'Invalid transport 2sigma band must be written as NaN')

    call sanitize_transport_output_triplet(.true., 3.0_dp, huge(1.0_dp), sanitized_mean, sanitized_std, sanitized_2sigma)

    call assert_close(sanitized_mean, 3.0_dp, 1.0d-12, 'Valid transport mean must stay unchanged')
    call assert_true(sanitized_std > 0.0_dp, 'Large valid transport standard deviation must remain positive after clipping')
    call assert_true(sanitized_2sigma > sanitized_std, 'Transport 2sigma band must exceed 1sigma after clipping')

end subroutine test_transport_output_sanitization_handles_invalid_sigma

subroutine test_local_reference_boundaries_avoid_axis_and_edge()

    use global_transport_fit_gorilla_mod, only: select_local_reference_boundaries

    integer, allocatable :: boundary_indices(:)
    integer :: i
    real(dp), allocatable :: boundary_s(:)

    allocate(boundary_s(27))
    boundary_s = [(real(i - 1, dp) / 26.0_dp, i=1, 27)]

    call select_local_reference_boundaries(boundary_s, 2, boundary_indices)

    call assert_true(all(boundary_indices >= 3), 'Local reference boundaries must stay away from the magnetic axis')
    call assert_true(all(boundary_indices <= 25), 'Local reference boundaries must stay away from the absorbing edge')

end subroutine test_local_reference_boundaries_avoid_axis_and_edge

subroutine test_local_flux_pair_recovery()

    real(dp) :: a_coeff
    real(dp) :: b_coeff
    real(dp) :: density_1
    real(dp) :: density_2
    real(dp) :: flux_1
    real(dp) :: flux_2
    real(dp) :: gradient_1
    real(dp) :: gradient_2
    real(dp) :: true_a
    real(dp) :: true_b

    true_a = 1.5d4
    true_b = 2.5d3
    gradient_1 = -6.0d0
    gradient_2 = 3.0d0
    density_1 = 2.0d13
    density_2 = 3.0d13
    flux_1 = density_1 * (true_a - true_b * gradient_1)
    flux_2 = density_2 * (true_a - true_b * gradient_2)

    call recover_transport_coefficients_from_flux_pair(flux_1, flux_2, density_1, density_2, gradient_1, gradient_2, &
        a_coeff, b_coeff)

    call assert_close(a_coeff, true_a, 1.0d-10, 'Local flux-pair A recovery failed')
    call assert_close(b_coeff, true_b, 1.0d-10, 'Local flux-pair B recovery failed')

end subroutine test_local_flux_pair_recovery

subroutine test_weighted_linear_transport_response()

    real(dp) :: a_coeff
    real(dp) :: a_std
    real(dp) :: b_coeff
    real(dp) :: b_std
    real(dp), allocatable :: gradients(:)
    real(dp), allocatable :: normalized_flux(:)
    real(dp), allocatable :: normalized_flux_variance(:)
    logical :: valid_fit

    gradients = (/-6.0_dp, -2.0_dp, 0.0_dp, 3.0_dp, 6.0_dp/)
    normalized_flux = 1.5d4 - 2.5d3 * gradients
    normalized_flux_variance = (/1.0d-2, 2.0d-2, 1.5d-2, 2.5d-2, 3.0d-2/)

    call fit_linear_transport_response(gradients, normalized_flux, normalized_flux_variance, a_coeff, b_coeff, a_std, b_std, &
        valid_fit)

    call assert_true(valid_fit, 'Weighted local transport regression must be valid')
    call assert_close(a_coeff, 1.5d4, 1.0d-10, 'Weighted local regression A recovery failed')
    call assert_close(b_coeff, 2.5d3, 1.0d-10, 'Weighted local regression B recovery failed')
    call assert_true(a_std > 0.0_dp, 'Weighted local regression A standard deviation must be positive')
    call assert_true(b_std > 0.0_dp, 'Weighted local regression B standard deviation must be positive')

end subroutine test_weighted_linear_transport_response

subroutine test_batch_mean_and_variance()

    real(dp), allocatable :: mean_values(:)
    real(dp), allocatable :: samples(:, :)
    real(dp), allocatable :: variance_of_mean(:)

    allocate(samples(2, 3))
    samples(:, 1) = (/1.0_dp, 4.0_dp/)
    samples(:, 2) = (/2.0_dp, 7.0_dp/)
    samples(:, 3) = (/3.0_dp, 10.0_dp/)

    call compute_mean_and_variance(samples, mean_values, variance_of_mean)

    call assert_close(mean_values(1), 2.0_dp, 1.0d-12, 'Batch mean for first observable is wrong')
    call assert_close(mean_values(2), 7.0_dp, 1.0d-12, 'Batch mean for second observable is wrong')
    call assert_close(variance_of_mean(1), 1.0_dp / 3.0_dp, 1.0d-12, 'Variance of mean for first observable is wrong')
    call assert_close(variance_of_mean(2), 3.0_dp, 1.0d-12, 'Variance of mean for second observable is wrong')

end subroutine test_batch_mean_and_variance

subroutine test_convergence_history_populated()

    type(global_transport_experiment_t), allocatable :: experiments(:)
    type(global_transport_fit_control_t) :: control
    type(global_transport_fit_result_t) :: result
    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: shell_volumes(:)
    real(dp), allocatable :: source_1(:)
    real(dp), allocatable :: source_2(:)
    real(dp), allocatable :: knot_s(:)
    real(dp), allocatable :: true_a(:)
    real(dp), allocatable :: true_b(:)

    call make_test_geometry(boundary_s, shell_volumes, source_1, source_2)

    control%n_knots_a = 3
    control%n_knots_b = 3
    control%max_lm_iterations = 25
    control%regularization_a = 1.0d-12
    control%regularization_b = 1.0d-12
    control%lm_damping = 1.0d-6

    call make_knot_points(boundary_s, control%n_knots_a, knot_s)
    allocate(true_a(control%n_knots_a))
    allocate(true_b(control%n_knots_b))
    true_a = (/5.0d-2, 2.0d-2, -1.0d-2/)
    true_b = log((/2.0d-2, 3.0d-2, 5.0d-2/))

    allocate(experiments(2))
    call generate_synthetic_experiment(boundary_s, shell_volumes, source_1, knot_s, knot_s, true_a, true_b, 0.0_dp, experiments(1))
    call generate_synthetic_experiment(boundary_s, shell_volumes, source_2, knot_s, knot_s, true_a, true_b, 0.0_dp, experiments(2))

    call fit_global_transport(experiments, control, result)

    call assert_true(allocated(result%history_objective), 'Convergence history not allocated')
    call assert_true(size(result%history_objective) > 0, 'Convergence history is empty')
    call assert_true(size(result%history_objective) == result%n_iterations, 'History length must match n_iterations')
    call assert_true(all(result%history_objective >= 0.0_dp), 'Objective must be non-negative')
    call assert_true(all(result%history_gradient_norm >= 0.0_dp), 'Gradient norm must be non-negative')
    call assert_true(all(result%history_damping > 0.0_dp), 'Damping must be positive')
    call assert_true(result%history_objective(size(result%history_objective)) < result%history_objective(1), &
        'Final objective must be smaller than initial objective')

end subroutine test_convergence_history_populated

subroutine compute_objective_only(experiments, control, basis, parameters, objective)

    type(global_transport_experiment_t), intent(in) :: experiments(:)
    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(in) :: basis(:, :)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(out) :: objective

    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: jacobian(:, :)

    call compute_objective_gradient_and_jacobian(experiments, control, basis, basis, parameters, objective, gradient, jacobian)

end subroutine compute_objective_only

subroutine make_test_geometry(boundary_s, shell_volumes, source_1, source_2)

    real(dp), allocatable, intent(out) :: boundary_s(:)
    real(dp), allocatable, intent(out) :: shell_volumes(:)
    real(dp), allocatable, intent(out) :: source_1(:)
    real(dp), allocatable, intent(out) :: source_2(:)

    boundary_s = (/0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp, 1.0_dp/)
    shell_volumes = (/1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp, 1.4_dp/)
    source_1 = (/1.5_dp, 0.8_dp, 0.2_dp, 0.0_dp, 0.0_dp/)
    source_2 = (/0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.6_dp/)

end subroutine make_test_geometry

subroutine make_knot_points(boundary_s, n_knots, knot_s)

    real(dp), intent(in) :: boundary_s(:)
    integer, intent(in) :: n_knots
    real(dp), allocatable, intent(out) :: knot_s(:)

    integer :: i
    real(dp) :: s_min
    real(dp) :: s_max

    allocate(knot_s(n_knots))
    s_min = boundary_s(1)
    s_max = boundary_s(size(boundary_s))
    do i = 1, n_knots
        knot_s(i) = s_min + real(i - 1, dp) * (s_max - s_min) / real(n_knots - 1, dp)
    end do

end subroutine make_knot_points

subroutine assert_close(value, expected, tolerance, message)

    real(dp), intent(in) :: value
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (abs(value - expected) > tolerance) then
        print *, trim(message), value, expected, tolerance
        stop 1
    end if

end subroutine assert_close

subroutine assert_true(condition, message)

    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
        print *, trim(message)
        stop 1
    end if

end subroutine assert_true

subroutine assert_relative_close(error_value, reference_value, tolerance, message)

    real(dp), intent(in) :: error_value
    real(dp), intent(in) :: reference_value
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (error_value / max(reference_value, 1.0d-12) > tolerance) then
        print *, trim(message), error_value, reference_value, tolerance
        stop 1
    end if

end subroutine assert_relative_close

end program test_global_transport_fit
