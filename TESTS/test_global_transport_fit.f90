program test_global_transport_fit

    use, intrinsic :: iso_fortran_env, only: dp => real64

    use global_transport_fit_core_mod, only: compute_objective_gradient_and_jacobian, fit_global_transport, &
        generate_synthetic_experiment
    use global_transport_fit_math_mod, only: build_piecewise_linear_basis
    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_fit_control_t, &
        global_transport_fit_result_t

    implicit none

    call test_fit_recovers_manufactured_profiles()
    call test_gradient_matches_finite_difference()

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
