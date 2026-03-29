module global_transport_fit_core_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    use global_transport_fit_math_mod, only: build_piecewise_linear_basis, build_second_difference_matrix, &
        compute_geometry_from_boundaries, solve_linear_system
    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_fit_control_t, &
        global_transport_fit_result_t
    use transport_statistics_mod, only: build_density_support_mask, build_inverse_variance_weights, &
        build_supported_boundary_mask

    implicit none

    private

    public :: compute_objective_gradient_and_jacobian
    public :: fit_global_transport
    public :: generate_synthetic_experiment
    public :: solve_forward_problem

contains

subroutine sanitize_parameters(control, parameters_in, parameters_out)

    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(in) :: parameters_in(:)
    real(dp), allocatable, intent(out) :: parameters_out(:)

    real(dp), parameter :: max_abs_log_b = 20.0_dp

    allocate(parameters_out(size(parameters_in)))
    parameters_out = parameters_in
    parameters_out(control%n_knots_a + 1:control%n_knots_a + control%n_knots_b) = max( &
        min(parameters_out(control%n_knots_a + 1:control%n_knots_a + control%n_knots_b), max_abs_log_b), -max_abs_log_b)

end subroutine sanitize_parameters

subroutine fit_global_transport(experiments, control, result)

    type(global_transport_experiment_t), intent(in) :: experiments(:)
    type(global_transport_fit_control_t), intent(in) :: control
    type(global_transport_fit_result_t), intent(out) :: result

    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: knot_s_a(:)
    real(dp), allocatable :: knot_s_b(:)
    real(dp), allocatable :: basis_a(:, :)
    real(dp), allocatable :: basis_b(:, :)
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: sanitized_parameters(:)
    real(dp), allocatable :: step(:)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: jacobian(:, :)
    real(dp), allocatable :: hessian(:, :)
    real(dp), allocatable :: covariance(:, :)
    logical :: converged
    real(dp) :: damping
    real(dp) :: gradient_norm
    real(dp) :: objective
    real(dp) :: step_norm
    real(dp) :: trial_objective
    integer :: n_parameters
    integer :: iteration
    logical :: accepted

    boundary_s = experiments(1)%boundary_s
    call build_knot_positions(boundary_s, control%n_knots_a, knot_s_a)
    call build_knot_positions(boundary_s, control%n_knots_b, knot_s_b)
    call build_piecewise_linear_basis(boundary_s, knot_s_a, basis_a)
    call build_piecewise_linear_basis(boundary_s, knot_s_b, basis_b)

    n_parameters = control%n_knots_a + control%n_knots_b
    allocate(parameters(n_parameters))
    call initialize_parameters(experiments, control, parameters)
    damping = control%lm_damping
    converged = .false.

    allocate(result%history_objective(control%max_lm_iterations))
    allocate(result%history_gradient_norm(control%max_lm_iterations))
    allocate(result%history_step_norm(control%max_lm_iterations))
    allocate(result%history_damping(control%max_lm_iterations))
    allocate(result%history_accepted(control%max_lm_iterations))

    call sanitize_parameters(control, parameters, sanitized_parameters)
    call move_alloc(sanitized_parameters, parameters)

    do iteration = 1, control%max_lm_iterations
        call compute_objective_gradient_and_jacobian(experiments, control, basis_a, basis_b, parameters, objective, gradient, jacobian)
        gradient_norm = maxval(abs(gradient))

        result%history_objective(iteration) = objective
        result%history_gradient_norm(iteration) = gradient_norm
        result%history_damping(iteration) = damping
        result%history_step_norm(iteration) = 0.0_dp
        result%history_accepted(iteration) = .false.

        if (gradient_norm < control%gradient_tolerance) then
            converged = .true.
            exit
        end if

        call build_lm_hessian(control, jacobian, parameters, hessian)
        call add_diagonal_shift(hessian, damping)
        call solve_step(hessian, gradient, step)

        step_norm = maxval(abs(step))
        result%history_step_norm(iteration) = step_norm

        if (step_norm < control%step_tolerance) then
            converged = .true.
            exit
        end if

        call try_lm_step(experiments, control, basis_a, basis_b, parameters, step, objective, damping, trial_objective, accepted)
        result%history_accepted(iteration) = accepted
        if (accepted) then
            call sanitize_parameters(control, parameters + step, sanitized_parameters)
            call move_alloc(sanitized_parameters, parameters)
            damping = max(control%lm_damping * 1.0d-12, damping * control%lm_damping_decrease)
            if (step_norm < 10.0_dp * control%step_tolerance) then
                converged = .true.
                exit
            end if
        else
            damping = damping * control%lm_damping_increase
        end if
    end do

    result%history_objective = result%history_objective(:min(iteration, control%max_lm_iterations))
    result%history_gradient_norm = result%history_gradient_norm(:min(iteration, control%max_lm_iterations))
    result%history_step_norm = result%history_step_norm(:min(iteration, control%max_lm_iterations))
    result%history_damping = result%history_damping(:min(iteration, control%max_lm_iterations))
    result%history_accepted = result%history_accepted(:min(iteration, control%max_lm_iterations))

    call compute_objective_gradient_and_jacobian(experiments, control, basis_a, basis_b, parameters, objective, gradient, jacobian)
    call build_lm_hessian(control, jacobian, parameters, hessian)
    call solve_covariance(hessian, covariance)
    call fill_result(boundary_s, knot_s_a, knot_s_b, basis_a, basis_b, parameters, gradient, covariance, objective, iteration, &
        converged, control, result)

end subroutine fit_global_transport

subroutine initialize_parameters(experiments, control, parameters)

    type(global_transport_experiment_t), intent(in) :: experiments(:)
    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(out) :: parameters(:)

    real(dp), allocatable :: cell_centers(:)
    real(dp), allocatable :: boundary_areas(:)
    real(dp), allocatable :: density_gradient(:)
    real(dp) :: average_flux
    real(dp) :: average_gradient
    real(dp) :: b_scale
    integer :: i

    parameters = 0.0_dp
    call compute_geometry_from_boundaries(experiments(1)%boundary_s, experiments(1)%shell_volumes, cell_centers, boundary_areas)
    allocate(density_gradient(size(experiments(1)%density) - 1))
    do i = 1, size(density_gradient)
        density_gradient(i) = (experiments(1)%density(i + 1) - experiments(1)%density(i)) / &
            max(cell_centers(i + 1) - cell_centers(i), 1.0d-12)
    end do
    average_flux = sum(abs(experiments(1)%flux(2:size(experiments(1)%flux))))
    average_flux = average_flux / real(max(size(experiments(1)%flux) - 1, 1), dp)
    average_gradient = sum(abs(density_gradient)) / real(max(size(density_gradient), 1), dp)
    b_scale = max(1.0d-4, average_flux / max(average_gradient, 1.0d-8))
    parameters(control%n_knots_a + 1:) = log(b_scale)

end subroutine initialize_parameters

subroutine generate_synthetic_experiment(boundary_s, shell_volumes, source, knot_s_a, knot_s_b, a_coefficients, b_coefficients, &
    outer_boundary_density, experiment)

    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: shell_volumes(:)
    real(dp), intent(in) :: source(:)
    real(dp), intent(in) :: knot_s_a(:)
    real(dp), intent(in) :: knot_s_b(:)
    real(dp), intent(in) :: a_coefficients(:)
    real(dp), intent(in) :: b_coefficients(:)
    real(dp), intent(in) :: outer_boundary_density
    type(global_transport_experiment_t), intent(out) :: experiment

    real(dp), allocatable :: basis_a(:, :)
    real(dp), allocatable :: basis_b(:, :)
    real(dp), allocatable :: a_boundary(:)
    real(dp), allocatable :: b_boundary(:)
    real(dp), allocatable :: boundary_areas(:)
    real(dp), allocatable :: cell_centers(:)
    real(dp), allocatable :: density(:)
    real(dp), allocatable :: flux(:)
    real(dp), allocatable :: integrated_flux(:)

    call build_piecewise_linear_basis(boundary_s, knot_s_a, basis_a)
    call build_piecewise_linear_basis(boundary_s, knot_s_b, basis_b)
    allocate(a_boundary(size(boundary_s)))
    allocate(b_boundary(size(boundary_s)))
    a_boundary = matmul(basis_a, a_coefficients)
    b_boundary = exp(matmul(basis_b, b_coefficients))

    call compute_geometry_from_boundaries(boundary_s, shell_volumes, cell_centers, boundary_areas)
    call solve_forward_problem(boundary_s, cell_centers, boundary_areas, shell_volumes, source, a_boundary, b_boundary, &
        outer_boundary_density, density, flux)

    experiment%boundary_s = boundary_s
    experiment%shell_volumes = shell_volumes
    experiment%source = source
    experiment%density = density
    experiment%flux = flux
    integrated_flux = boundary_areas * flux
    experiment%integrated_flux = integrated_flux
    allocate(experiment%density_variance(size(density)))
    allocate(experiment%flux_variance(size(flux)))
    allocate(experiment%integrated_flux_variance(size(flux)))
    experiment%density_variance = 1.0d-10
    experiment%flux_variance = 1.0d-10
    experiment%integrated_flux_variance = max(boundary_areas**2 * experiment%flux_variance, 1.0d-10)

end subroutine generate_synthetic_experiment

subroutine build_knot_positions(boundary_s, n_knots, knot_s)

    real(dp), intent(in) :: boundary_s(:)
    integer, intent(in) :: n_knots
    real(dp), allocatable, intent(out) :: knot_s(:)

    integer :: i
    real(dp) :: s_max
    real(dp) :: s_min

    allocate(knot_s(n_knots))
    s_min = boundary_s(1)
    s_max = boundary_s(size(boundary_s))
    if (n_knots == 1) then
        knot_s(1) = 0.5_dp * (s_min + s_max)
        return
    end if

    do i = 1, n_knots
        knot_s(i) = s_min + real(i - 1, dp) * (s_max - s_min) / real(n_knots - 1, dp)
    end do

end subroutine build_knot_positions

subroutine compute_objective_gradient_and_jacobian(experiments, control, basis_a, basis_b, parameters, objective, gradient, jacobian)

    type(global_transport_experiment_t), intent(in) :: experiments(:)
    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(in) :: basis_a(:, :)
    real(dp), intent(in) :: basis_b(:, :)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(out) :: objective
    real(dp), allocatable, intent(out) :: gradient(:)
    real(dp), allocatable, intent(out) :: jacobian(:, :)

    real(dp), allocatable :: regularization_matrix_a(:, :)
    real(dp), allocatable :: regularization_matrix_b(:, :)
    real(dp), allocatable :: weighted_jacobian(:, :)
    integer :: row_offset

    call build_second_difference_matrix(control%n_knots_a, regularization_matrix_a)
    call build_second_difference_matrix(control%n_knots_b, regularization_matrix_b)

    allocate(gradient(size(parameters)))
    allocate(jacobian(total_observable_count(experiments, control%use_flux_objective), size(parameters)))
    jacobian = 0.0_dp
    gradient = 0.0_dp
    objective = 0.0_dp
    row_offset = 0

    call accumulate_experiment_terms(experiments, control, basis_a, basis_b, parameters, objective, gradient, jacobian, row_offset)
    call add_regularization_terms(control, regularization_matrix_a, regularization_matrix_b, parameters, objective, gradient)
    call apply_weighting_to_jacobian(experiments, control, jacobian, weighted_jacobian)
    jacobian = weighted_jacobian

end subroutine compute_objective_gradient_and_jacobian

subroutine accumulate_experiment_terms(experiments, control, basis_a, basis_b, parameters, objective, gradient, jacobian, row_offset)

    type(global_transport_experiment_t), intent(in) :: experiments(:)
    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(in) :: basis_a(:, :)
    real(dp), intent(in) :: basis_b(:, :)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(inout) :: objective
    real(dp), intent(inout) :: gradient(:)
    real(dp), intent(inout) :: jacobian(:, :)
    integer, intent(inout) :: row_offset

    integer :: i

    do i = 1, size(experiments)
        call accumulate_single_experiment(experiments(i), control, basis_a, basis_b, parameters, objective, gradient, jacobian, &
            row_offset)
    end do

end subroutine accumulate_experiment_terms

subroutine accumulate_single_experiment(experiment, control, basis_a, basis_b, parameters, objective, gradient, jacobian, row_offset)

    type(global_transport_experiment_t), intent(in) :: experiment
    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(in) :: basis_a(:, :)
    real(dp), intent(in) :: basis_b(:, :)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(inout) :: objective
    real(dp), intent(inout) :: gradient(:)
    real(dp), intent(inout) :: jacobian(:, :)
    integer, intent(inout) :: row_offset

    real(dp), allocatable :: a_boundary(:)
    real(dp), allocatable :: b_boundary(:)
    real(dp), allocatable :: da_dp(:, :)
    real(dp), allocatable :: db_dp(:, :)
    real(dp), allocatable :: cell_centers(:)
    real(dp), allocatable :: boundary_areas(:)
    real(dp), allocatable :: density(:)
    real(dp), allocatable :: flux(:)
    real(dp), allocatable :: integrated_flux(:)
    real(dp), allocatable :: system_matrix(:, :)
    real(dp), allocatable :: flux_n(:, :)
    real(dp), allocatable :: flux_param_partial(:, :)
    real(dp), allocatable :: integrated_flux_n(:, :)
    real(dp), allocatable :: integrated_flux_param_partial(:, :)
    real(dp), allocatable :: density_sensitivity(:, :)
    real(dp), allocatable :: observable_jacobian(:, :)
    real(dp), allocatable :: density_weights(:)
    real(dp), allocatable :: flux_weights(:)
    real(dp), allocatable :: residual_density(:)
    real(dp), allocatable :: residual_flux(:)
    real(dp), allocatable :: adjoint_rhs(:)
    real(dp), allocatable :: adjoint_state(:)
    real(dp), allocatable :: divergence_matrix(:, :)
    real(dp) :: experiment_objective

    call evaluate_profiles(parameters, basis_a, basis_b, a_boundary, b_boundary, da_dp, db_dp)
    call compute_geometry_from_boundaries(experiment%boundary_s, experiment%shell_volumes, cell_centers, boundary_areas)
    call solve_forward_with_partials(experiment, control%outer_boundary_density, a_boundary, b_boundary, cell_centers, boundary_areas, &
        density, flux, system_matrix, flux_n, flux_param_partial, divergence_matrix, da_dp, db_dp)
    call multiply_flux_by_geometry(boundary_areas, flux_n, integrated_flux_n)
    call multiply_flux_param_by_geometry(boundary_areas, flux_param_partial, integrated_flux_param_partial)
    integrated_flux = boundary_areas * flux
    call build_experiment_weights(control, experiment, density_weights, flux_weights)

    residual_density = density - experiment%density
    residual_flux = integrated_flux - experiment%integrated_flux
    experiment_objective = 0.5_dp * sum(density_weights * residual_density**2)
    if (control%use_flux_objective) then
        experiment_objective = experiment_objective + 0.5_dp * sum(flux_weights * residual_flux**2)
    end if
    objective = objective + experiment_objective

    call build_adjoint_rhs(integrated_flux_n, density_weights, flux_weights, residual_density, residual_flux, &
        control%use_flux_objective, adjoint_rhs)
    call solve_adjoint(system_matrix, adjoint_rhs, adjoint_state)
    gradient = gradient - matmul(transpose(matmul(divergence_matrix, flux_param_partial)), adjoint_state)
    if (control%use_flux_objective) then
        gradient = gradient + matmul(transpose(integrated_flux_param_partial), flux_weights * residual_flux)
    end if

    call build_density_sensitivity(system_matrix, divergence_matrix, flux_param_partial, density_sensitivity)
    call build_observable_jacobian(integrated_flux_n, integrated_flux_param_partial, density_sensitivity, control%use_flux_objective, &
        observable_jacobian)
    jacobian(row_offset + 1:row_offset + size(observable_jacobian, 1), :) = observable_jacobian
    row_offset = row_offset + size(observable_jacobian, 1)

end subroutine accumulate_single_experiment

subroutine evaluate_profiles(parameters, basis_a, basis_b, a_boundary, b_boundary, da_dp, db_dp)

    real(dp), intent(in) :: parameters(:)
    real(dp), intent(in) :: basis_a(:, :)
    real(dp), intent(in) :: basis_b(:, :)
    real(dp), allocatable, intent(out) :: a_boundary(:)
    real(dp), allocatable, intent(out) :: b_boundary(:)
    real(dp), allocatable, intent(out) :: da_dp(:, :)
    real(dp), allocatable, intent(out) :: db_dp(:, :)

    integer :: n_a
    integer :: n_b
    integer :: n_parameters
    real(dp), allocatable :: log_b(:)

    n_a = size(basis_a, 2)
    n_b = size(basis_b, 2)
    n_parameters = size(parameters)
    allocate(a_boundary(size(basis_a, 1)))
    allocate(b_boundary(size(basis_b, 1)))
    allocate(da_dp(size(a_boundary), n_parameters))
    allocate(db_dp(size(b_boundary), n_parameters))
    allocate(log_b(size(b_boundary)))

    a_boundary = matmul(basis_a, parameters(1:n_a))
    log_b = max(min(matmul(basis_b, parameters(n_a + 1:n_a + n_b)), 20.0_dp), -20.0_dp)
    b_boundary = exp(log_b)

    da_dp = 0.0_dp
    da_dp(:, 1:n_a) = basis_a
    db_dp = 0.0_dp
    db_dp(:, n_a + 1:n_a + n_b) = spread(b_boundary, 2, n_b) * basis_b

end subroutine evaluate_profiles

subroutine solve_forward_with_partials(experiment, outer_boundary_density, a_boundary, b_boundary, cell_centers, boundary_areas, &
    density, flux, system_matrix, flux_n, flux_param_partial, divergence_matrix, da_dp, db_dp)

    type(global_transport_experiment_t), intent(in) :: experiment
    real(dp), intent(in) :: outer_boundary_density
    real(dp), intent(in) :: a_boundary(:)
    real(dp), intent(in) :: b_boundary(:)
    real(dp), intent(in) :: cell_centers(:)
    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: da_dp(:, :)
    real(dp), intent(in) :: db_dp(:, :)
    real(dp), allocatable, intent(out) :: density(:)
    real(dp), allocatable, intent(out) :: flux(:)
    real(dp), allocatable, intent(out) :: system_matrix(:, :)
    real(dp), allocatable, intent(out) :: flux_n(:, :)
    real(dp), allocatable, intent(out) :: flux_param_partial(:, :)
    real(dp), allocatable, intent(out) :: divergence_matrix(:, :)

    real(dp), allocatable :: flux_const(:)
    real(dp), allocatable :: rhs(:, :)
    real(dp), allocatable :: density_matrix(:, :)
    call build_flux_operators(experiment%boundary_s, cell_centers, a_boundary, b_boundary, outer_boundary_density, flux_n, flux_const)
    call build_divergence_matrix(boundary_areas, divergence_matrix)
    system_matrix = matmul(divergence_matrix, flux_n)
    allocate(rhs(size(experiment%source), 1))
    rhs(:, 1) = experiment%shell_volumes * experiment%source - matmul(divergence_matrix, flux_const)
    call solve_linear_system(system_matrix, rhs, density_matrix)
    allocate(density(size(experiment%source)))
    density = density_matrix(:, 1)

    flux = matmul(flux_n, density) + flux_const
    call build_flux_param_partial(experiment%boundary_s, cell_centers, density, outer_boundary_density, da_dp, db_dp, &
        flux_param_partial)

end subroutine solve_forward_with_partials

subroutine solve_forward_problem(boundary_s, cell_centers, boundary_areas, shell_volumes, source, a_boundary, b_boundary, &
    outer_boundary_density, density, flux)

    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: cell_centers(:)
    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: shell_volumes(:)
    real(dp), intent(in) :: source(:)
    real(dp), intent(in) :: a_boundary(:)
    real(dp), intent(in) :: b_boundary(:)
    real(dp), intent(in) :: outer_boundary_density
    real(dp), allocatable, intent(out) :: density(:)
    real(dp), allocatable, intent(out) :: flux(:)

    type(global_transport_experiment_t) :: experiment
    real(dp), allocatable :: flux_n(:, :)
    real(dp), allocatable :: flux_param_partial(:, :)
    real(dp), allocatable :: divergence_matrix(:, :)
    real(dp), allocatable :: system_matrix(:, :)
    real(dp), allocatable :: da_dp(:, :)
    real(dp), allocatable :: db_dp(:, :)

    experiment%boundary_s = boundary_s
    experiment%shell_volumes = shell_volumes
    experiment%source = source
    allocate(da_dp(size(boundary_s), 0))
    allocate(db_dp(size(boundary_s), 0))
    call solve_forward_with_partials(experiment, outer_boundary_density, a_boundary, b_boundary, cell_centers, boundary_areas, &
        density, flux, system_matrix, flux_n, flux_param_partial, divergence_matrix, da_dp, db_dp)

end subroutine solve_forward_problem

subroutine build_flux_operators(boundary_s, cell_centers, a_boundary, b_boundary, outer_boundary_density, flux_n, flux_const)

    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: cell_centers(:)
    real(dp), intent(in) :: a_boundary(:)
    real(dp), intent(in) :: b_boundary(:)
    real(dp), intent(in) :: outer_boundary_density
    real(dp), allocatable, intent(out) :: flux_n(:, :)
    real(dp), allocatable, intent(out) :: flux_const(:)

    integer :: j
    integer :: n_boundaries
    integer :: n_cells
    real(dp) :: delta_s

    n_boundaries = size(boundary_s)
    n_cells = n_boundaries - 1
    allocate(flux_n(n_boundaries, n_cells))
    allocate(flux_const(n_boundaries))
    flux_n = 0.0_dp
    flux_const = 0.0_dp

    do j = 2, n_cells
        delta_s = cell_centers(j) - cell_centers(j - 1)
        flux_n(j, j - 1) = 0.5_dp * a_boundary(j) + b_boundary(j) / delta_s
        flux_n(j, j) = 0.5_dp * a_boundary(j) - b_boundary(j) / delta_s
    end do

    delta_s = boundary_s(n_boundaries) - cell_centers(n_cells)
    flux_n(n_boundaries, n_cells) = 0.5_dp * a_boundary(n_boundaries) + b_boundary(n_boundaries) / delta_s
    flux_const(n_boundaries) = (0.5_dp * a_boundary(n_boundaries) - b_boundary(n_boundaries) / delta_s) * outer_boundary_density

end subroutine build_flux_operators

subroutine build_flux_param_partial(boundary_s, cell_centers, density, outer_boundary_density, da_dp, db_dp, flux_param_partial)

    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: cell_centers(:)
    real(dp), intent(in) :: density(:)
    real(dp), intent(in) :: outer_boundary_density
    real(dp), intent(in) :: da_dp(:, :)
    real(dp), intent(in) :: db_dp(:, :)
    real(dp), allocatable, intent(out) :: flux_param_partial(:, :)

    integer :: j
    integer :: n_boundaries
    integer :: n_cells
    real(dp) :: delta_s
    real(dp) :: dflux_da
    real(dp) :: dflux_db

    n_boundaries = size(boundary_s)
    n_cells = size(density)
    allocate(flux_param_partial(n_boundaries, size(da_dp, 2)))
    flux_param_partial = 0.0_dp

    do j = 2, n_cells
        delta_s = cell_centers(j) - cell_centers(j - 1)
        dflux_da = 0.5_dp * (density(j - 1) + density(j))
        dflux_db = -(density(j) - density(j - 1)) / delta_s
        flux_param_partial(j, :) = dflux_da * da_dp(j, :) + dflux_db * db_dp(j, :)
    end do

    delta_s = boundary_s(n_boundaries) - cell_centers(n_cells)
    dflux_da = 0.5_dp * (density(n_cells) + outer_boundary_density)
    dflux_db = -(outer_boundary_density - density(n_cells)) / delta_s
    flux_param_partial(n_boundaries, :) = dflux_da * da_dp(n_boundaries, :) + dflux_db * db_dp(n_boundaries, :)

end subroutine build_flux_param_partial

subroutine build_divergence_matrix(boundary_areas, divergence_matrix)

    real(dp), intent(in) :: boundary_areas(:)
    real(dp), allocatable, intent(out) :: divergence_matrix(:, :)

    integer :: i
    integer :: n_cells

    n_cells = size(boundary_areas) - 1
    allocate(divergence_matrix(n_cells, n_cells + 1))
    divergence_matrix = 0.0_dp

    do i = 1, n_cells
        divergence_matrix(i, i) = -boundary_areas(i)
        divergence_matrix(i, i + 1) = boundary_areas(i + 1)
    end do

end subroutine build_divergence_matrix

subroutine multiply_flux_by_geometry(boundary_areas, flux_operator, integrated_flux_operator)

    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: flux_operator(:, :)
    real(dp), allocatable, intent(out) :: integrated_flux_operator(:, :)

    integer :: i

    allocate(integrated_flux_operator(size(flux_operator, 1), size(flux_operator, 2)))
    do i = 1, size(boundary_areas)
        integrated_flux_operator(i, :) = boundary_areas(i) * flux_operator(i, :)
    end do

end subroutine multiply_flux_by_geometry

subroutine multiply_flux_param_by_geometry(boundary_areas, flux_param_partial, integrated_flux_param_partial)

    real(dp), intent(in) :: boundary_areas(:)
    real(dp), intent(in) :: flux_param_partial(:, :)
    real(dp), allocatable, intent(out) :: integrated_flux_param_partial(:, :)

    integer :: i

    allocate(integrated_flux_param_partial(size(flux_param_partial, 1), size(flux_param_partial, 2)))
    do i = 1, size(boundary_areas)
        integrated_flux_param_partial(i, :) = boundary_areas(i) * flux_param_partial(i, :)
    end do

end subroutine multiply_flux_param_by_geometry

subroutine build_experiment_weights(control, experiment, density_weights, flux_weights)

    type(global_transport_fit_control_t), intent(in) :: control
    type(global_transport_experiment_t), intent(in) :: experiment
    real(dp), allocatable, intent(out) :: density_weights(:)
    real(dp), allocatable, intent(out) :: flux_weights(:)

    logical, allocatable :: density_supported(:)
    logical, allocatable :: flux_supported(:)

    call build_density_support_mask(experiment%source, experiment%density, control%density_support_fraction, density_supported)
    call build_supported_boundary_mask(experiment%integrated_flux, experiment%integrated_flux_variance, &
        control%support_sigma_multiplier, flux_supported)
    call build_inverse_variance_weights(experiment%density, experiment%density_variance, control%density_relative_sigma_floor, &
        density_supported, density_weights)
    call build_inverse_variance_weights(experiment%integrated_flux, experiment%integrated_flux_variance, &
        control%flux_relative_sigma_floor, flux_supported, flux_weights)

end subroutine build_experiment_weights

subroutine build_adjoint_rhs(flux_n, density_weights, flux_weights, residual_density, residual_flux, use_flux_objective, &
    adjoint_rhs)

    real(dp), intent(in) :: flux_n(:, :)
    real(dp), intent(in) :: density_weights(:)
    real(dp), intent(in) :: flux_weights(:)
    real(dp), intent(in) :: residual_density(:)
    real(dp), intent(in) :: residual_flux(:)
    logical, intent(in) :: use_flux_objective
    real(dp), allocatable, intent(out) :: adjoint_rhs(:)

    allocate(adjoint_rhs(size(residual_density)))
    adjoint_rhs = density_weights * residual_density
    if (use_flux_objective) then
        adjoint_rhs = adjoint_rhs + matmul(transpose(flux_n), flux_weights * residual_flux)
    end if

end subroutine build_adjoint_rhs

subroutine solve_adjoint(system_matrix, adjoint_rhs, adjoint_state)

    real(dp), intent(in) :: system_matrix(:, :)
    real(dp), intent(in) :: adjoint_rhs(:)
    real(dp), allocatable, intent(out) :: adjoint_state(:)

    real(dp), allocatable :: rhs(:, :)
    real(dp), allocatable :: solution(:, :)
    real(dp), allocatable :: system_transpose(:, :)

    allocate(rhs(size(adjoint_rhs), 1))
    allocate(system_transpose(size(system_matrix, 2), size(system_matrix, 1)))
    rhs(:, 1) = adjoint_rhs
    system_transpose = transpose(system_matrix)
    call solve_linear_system(system_transpose, rhs, solution)
    allocate(adjoint_state(size(adjoint_rhs)))
    adjoint_state = solution(:, 1)

end subroutine solve_adjoint

subroutine build_density_sensitivity(system_matrix, divergence_matrix, flux_param_partial, density_sensitivity)

    real(dp), intent(in) :: system_matrix(:, :)
    real(dp), intent(in) :: divergence_matrix(:, :)
    real(dp), intent(in) :: flux_param_partial(:, :)
    real(dp), allocatable, intent(out) :: density_sensitivity(:, :)

    real(dp), allocatable :: rhs(:, :)

    rhs = -matmul(divergence_matrix, flux_param_partial)
    call solve_linear_system(system_matrix, rhs, density_sensitivity)

end subroutine build_density_sensitivity

subroutine build_observable_jacobian(flux_n, flux_param_partial, density_sensitivity, use_flux_objective, observable_jacobian)

    real(dp), intent(in) :: flux_n(:, :)
    real(dp), intent(in) :: flux_param_partial(:, :)
    real(dp), intent(in) :: density_sensitivity(:, :)
    logical, intent(in) :: use_flux_objective
    real(dp), allocatable, intent(out) :: observable_jacobian(:, :)

    integer :: n_density
    integer :: n_flux

    n_density = size(density_sensitivity, 1)
    n_flux = 0
    if (use_flux_objective) n_flux = size(flux_param_partial, 1)
    allocate(observable_jacobian(n_density + n_flux, size(density_sensitivity, 2)))
    observable_jacobian(1:n_density, :) = density_sensitivity
    if (use_flux_objective) then
        observable_jacobian(n_density + 1:n_density + n_flux, :) = matmul(flux_n, density_sensitivity) + flux_param_partial
    end if

end subroutine build_observable_jacobian

subroutine add_regularization_terms(control, regularization_matrix_a, regularization_matrix_b, parameters, objective, gradient)

    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(in) :: regularization_matrix_a(:, :)
    real(dp), intent(in) :: regularization_matrix_b(:, :)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(inout) :: objective
    real(dp), intent(inout) :: gradient(:)

    real(dp), allocatable :: coeff_a(:)
    real(dp), allocatable :: coeff_b(:)
    real(dp), allocatable :: diff_a(:)
    real(dp), allocatable :: diff_b(:)

    coeff_a = parameters(1:control%n_knots_a)
    coeff_b = parameters(control%n_knots_a + 1:control%n_knots_a + control%n_knots_b)
    diff_a = matmul(regularization_matrix_a, coeff_a)
    diff_b = matmul(regularization_matrix_b, coeff_b)

    objective = objective + 0.5_dp * control%regularization_a * sum(diff_a**2) + &
        0.5_dp * control%regularization_b * sum(diff_b**2)
    if (size(regularization_matrix_a, 1) > 0) then
        gradient(1:control%n_knots_a) = gradient(1:control%n_knots_a) + &
            control%regularization_a * matmul(transpose(regularization_matrix_a), diff_a)
    end if
    if (size(regularization_matrix_b, 1) > 0) then
        gradient(control%n_knots_a + 1:control%n_knots_a + control%n_knots_b) = &
            gradient(control%n_knots_a + 1:control%n_knots_a + control%n_knots_b) + &
            control%regularization_b * matmul(transpose(regularization_matrix_b), diff_b)
    end if

end subroutine add_regularization_terms

subroutine apply_weighting_to_jacobian(experiments, control, jacobian, weighted_jacobian)

    type(global_transport_experiment_t), intent(in) :: experiments(:)
    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(in) :: jacobian(:, :)
    real(dp), allocatable, intent(out) :: weighted_jacobian(:, :)

    integer :: i
    integer :: n_density
    integer :: n_flux
    integer :: row_offset
    real(dp), allocatable :: density_weights(:)
    real(dp), allocatable :: flux_weights(:)

    allocate(weighted_jacobian(size(jacobian, 1), size(jacobian, 2)))
    row_offset = 0
    do i = 1, size(experiments)
        call build_experiment_weights(control, experiments(i), density_weights, flux_weights)
        n_density = size(experiments(i)%density)
        n_flux = 0
        if (control%use_flux_objective) n_flux = size(experiments(i)%integrated_flux)
        weighted_jacobian(row_offset + 1:row_offset + n_density, :) = &
            spread(sqrt(density_weights), 2, size(jacobian, 2)) * jacobian(row_offset + 1:row_offset + n_density, :)
        if (n_flux > 0) then
            weighted_jacobian(row_offset + n_density + 1:row_offset + n_density + n_flux, :) = &
                spread(sqrt(flux_weights), 2, size(jacobian, 2)) * &
                jacobian(row_offset + n_density + 1:row_offset + n_density + n_flux, :)
        end if
        row_offset = row_offset + n_density + n_flux
    end do

end subroutine apply_weighting_to_jacobian

subroutine build_lm_hessian(control, jacobian, parameters, hessian)

    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(in) :: jacobian(:, :)
    real(dp), intent(in) :: parameters(:)
    real(dp), allocatable, intent(out) :: hessian(:, :)

    real(dp), allocatable :: regularization_matrix_a(:, :)
    real(dp), allocatable :: regularization_matrix_b(:, :)

    hessian = matmul(transpose(jacobian), jacobian)
    call build_second_difference_matrix(control%n_knots_a, regularization_matrix_a)
    call build_second_difference_matrix(control%n_knots_b, regularization_matrix_b)
    if (size(regularization_matrix_a, 1) > 0) then
        hessian(1:control%n_knots_a, 1:control%n_knots_a) = hessian(1:control%n_knots_a, 1:control%n_knots_a) + &
            control%regularization_a * matmul(transpose(regularization_matrix_a), regularization_matrix_a)
    end if
    if (size(regularization_matrix_b, 1) > 0) then
        hessian(control%n_knots_a + 1:, control%n_knots_a + 1:) = hessian(control%n_knots_a + 1:, control%n_knots_a + 1:) + &
            control%regularization_b * matmul(transpose(regularization_matrix_b), regularization_matrix_b)
    end if

end subroutine build_lm_hessian

subroutine add_diagonal_shift(hessian, damping)

    real(dp), intent(inout) :: hessian(:, :)
    real(dp), intent(in) :: damping

    integer :: i

    do i = 1, size(hessian, 1)
        hessian(i, i) = hessian(i, i) + damping
    end do

end subroutine add_diagonal_shift

subroutine solve_step(hessian, gradient, step)

    real(dp), intent(in) :: hessian(:, :)
    real(dp), intent(in) :: gradient(:)
    real(dp), allocatable, intent(out) :: step(:)

    real(dp), allocatable :: rhs(:, :)
    real(dp), allocatable :: solution(:, :)

    allocate(rhs(size(gradient), 1))
    rhs(:, 1) = -gradient
    call solve_linear_system(hessian, rhs, solution)
    allocate(step(size(gradient)))
    step = solution(:, 1)

end subroutine solve_step

subroutine try_lm_step(experiments, control, basis_a, basis_b, parameters, step, objective, damping, trial_objective, accepted)

    type(global_transport_experiment_t), intent(in) :: experiments(:)
    type(global_transport_fit_control_t), intent(in) :: control
    real(dp), intent(in) :: basis_a(:, :)
    real(dp), intent(in) :: basis_b(:, :)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(in) :: step(:)
    real(dp), intent(in) :: objective
    real(dp), intent(in) :: damping
    real(dp), intent(out) :: trial_objective
    logical, intent(out) :: accepted

    real(dp), allocatable :: trial_parameters(:)
    real(dp), allocatable :: trial_gradient(:)
    real(dp), allocatable :: trial_jacobian(:, :)

    call sanitize_parameters(control, parameters + step, trial_parameters)
    call compute_objective_gradient_and_jacobian(experiments, control, basis_a, basis_b, trial_parameters, trial_objective, &
        trial_gradient, trial_jacobian)
    accepted = trial_objective < objective

end subroutine try_lm_step

subroutine solve_covariance(hessian, covariance)

    real(dp), intent(in) :: hessian(:, :)
    real(dp), allocatable, intent(out) :: covariance(:, :)

    real(dp), allocatable :: identity(:, :)
    integer :: i

    allocate(identity(size(hessian, 1), size(hessian, 1)))
    identity = 0.0_dp
    do i = 1, size(hessian, 1)
        identity(i, i) = 1.0_dp
    end do
    call solve_linear_system(hessian, identity, covariance)

end subroutine solve_covariance

subroutine fill_result(boundary_s, knot_s_a, knot_s_b, basis_a, basis_b, parameters, gradient, covariance, objective, iteration, &
    converged, control, result)

    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: knot_s_a(:)
    real(dp), intent(in) :: knot_s_b(:)
    real(dp), intent(in) :: basis_a(:, :)
    real(dp), intent(in) :: basis_b(:, :)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(in) :: gradient(:)
    real(dp), intent(in) :: covariance(:, :)
    real(dp), intent(in) :: objective
    integer, intent(in) :: iteration
    logical, intent(in) :: converged
    type(global_transport_fit_control_t), intent(in) :: control
    type(global_transport_fit_result_t), intent(inout) :: result

    real(dp), allocatable :: a_boundary(:)
    real(dp), allocatable :: b_boundary(:)
    real(dp), allocatable :: da_dp(:, :)
    real(dp), allocatable :: db_dp(:, :)
    integer :: i

    call evaluate_profiles(parameters, basis_a, basis_b, a_boundary, b_boundary, da_dp, db_dp)

    result%converged = converged .or. (maxval(abs(gradient)) < control%gradient_tolerance)
    result%n_iterations = iteration
    result%objective = objective
    result%parameters = parameters
    result%gradient = gradient
    result%covariance = covariance
    result%a_profile = a_boundary
    result%b_profile = b_boundary
    result%knot_s_a = knot_s_a
    result%knot_s_b = knot_s_b
    allocate(result%a_std(size(boundary_s)))
    allocate(result%b_std(size(boundary_s)))
    do i = 1, size(boundary_s)
        result%a_std(i) = sqrt(max(dot_product(da_dp(i, :), matmul(covariance, da_dp(i, :))), 0.0_dp))
        result%b_std(i) = sqrt(max(dot_product(db_dp(i, :), matmul(covariance, db_dp(i, :))), 0.0_dp))
    end do

end subroutine fill_result

integer function total_observable_count(experiments, use_flux_objective)

    type(global_transport_experiment_t), intent(in) :: experiments(:)
    logical, intent(in) :: use_flux_objective

    integer :: i

    total_observable_count = 0
    do i = 1, size(experiments)
        total_observable_count = total_observable_count + size(experiments(i)%density)
        if (use_flux_objective) total_observable_count = total_observable_count + size(experiments(i)%flux)
    end do

end function total_observable_count

end module global_transport_fit_core_mod
