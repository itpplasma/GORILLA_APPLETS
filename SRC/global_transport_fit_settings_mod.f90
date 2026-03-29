module global_transport_fit_settings_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    use global_transport_fit_types_mod, only: global_transport_fit_control_t

    implicit none

    private

    public :: load_global_transport_fit_inp

    type(global_transport_fit_control_t), public, protected :: control
    character(len=32), public, protected :: run_mode = 'files'
    character(len=256), public, protected :: filename_boundary_s = 'boundary_s.dat'
    character(len=256), public, protected :: filename_shell_volumes = 'shell_volumes.dat'
    character(len=256), public, protected :: filename_source_1 = 'source_1.dat'
    character(len=256), public, protected :: filename_density_1 = 'density_1.dat'
    character(len=256), public, protected :: filename_density_variance_1 = 'density_variance_1.dat'
    character(len=256), public, protected :: filename_flux_1 = 'flux_1.dat'
    character(len=256), public, protected :: filename_flux_variance_1 = 'flux_variance_1.dat'
    character(len=256), public, protected :: filename_source_2 = 'source_2.dat'
    character(len=256), public, protected :: filename_density_2 = 'density_2.dat'
    character(len=256), public, protected :: filename_density_variance_2 = 'density_variance_2.dat'
    character(len=256), public, protected :: filename_flux_2 = 'flux_2.dat'
    character(len=256), public, protected :: filename_flux_variance_2 = 'flux_variance_2.dat'
    character(len=256), public, protected :: filename_fit_summary = 'global_transport_fit_summary.dat'
    character(len=256), public, protected :: filename_fit_profiles = 'global_transport_fit_profiles.dat'
    character(len=256), public, protected :: filename_local_profiles = 'global_transport_local_profiles.dat'
    character(len=256), public, protected :: filename_comparison_profiles = 'global_transport_comparison_profiles.dat'
    character(len=256), public, protected :: filename_experiment_1_summary = 'global_transport_experiment_1.dat'
    character(len=256), public, protected :: filename_experiment_2_summary = 'global_transport_experiment_2.dat'
    integer, public, protected :: n_local_reference_surfaces = 5
    real(dp), public, protected :: local_density_gradient_1 = -6.0d0
    real(dp), public, protected :: local_density_gradient_2 = 6.0d0
    real(dp), public, protected :: local_density_profile_reference_s = 0.5d0
    real(dp), public, protected :: local_v_E = 0.0d0
    real(dp), public, protected :: source_gradient_1 = -6.0d0
    real(dp), public, protected :: source_gradient_2 = 6.0d0
    real(dp), public, protected :: source_reference_s_1 = 0.25d0
    real(dp), public, protected :: source_reference_s_2 = 0.75d0

contains

subroutine load_global_transport_fit_inp()

    integer :: fit_unit
    integer :: max_lm_iterations
    integer :: n_knots_a
    integer :: n_knots_b
    real(dp) :: gradient_tolerance
    real(dp) :: lm_damping
    real(dp) :: lm_damping_decrease
    real(dp) :: lm_damping_increase
    real(dp) :: outer_boundary_density
    real(dp) :: regularization_a
    real(dp) :: regularization_b
    real(dp) :: step_tolerance
    logical :: use_flux_objective

    namelist /global_transport_fit_nml/ run_mode, filename_boundary_s, filename_shell_volumes, filename_source_1, &
        filename_density_1, filename_density_variance_1, filename_flux_1, filename_flux_variance_1, &
        filename_source_2, filename_density_2, filename_density_variance_2, filename_flux_2, &
        filename_flux_variance_2, filename_fit_summary, filename_fit_profiles, filename_local_profiles, &
        filename_comparison_profiles, filename_experiment_1_summary, filename_experiment_2_summary, &
        n_local_reference_surfaces, local_density_gradient_1, local_density_gradient_2, local_density_profile_reference_s, &
        local_v_E, source_gradient_1, source_gradient_2, source_reference_s_1, source_reference_s_2, n_knots_a, &
        n_knots_b, use_flux_objective, regularization_a, regularization_b, max_lm_iterations, gradient_tolerance, &
        step_tolerance, lm_damping, lm_damping_increase, lm_damping_decrease, outer_boundary_density

    n_knots_a = control%n_knots_a
    n_knots_b = control%n_knots_b
    regularization_a = control%regularization_a
    regularization_b = control%regularization_b
    max_lm_iterations = control%max_lm_iterations
    gradient_tolerance = control%gradient_tolerance
    step_tolerance = control%step_tolerance
    lm_damping = control%lm_damping
    lm_damping_increase = control%lm_damping_increase
    lm_damping_decrease = control%lm_damping_decrease
    outer_boundary_density = control%outer_boundary_density
    use_flux_objective = control%use_flux_objective

    open(newunit=fit_unit, file='global_transport_fit.inp', status='old', action='read')
    read(fit_unit, nml=global_transport_fit_nml)
    close(fit_unit)

    control%n_knots_a = n_knots_a
    control%n_knots_b = n_knots_b
    control%regularization_a = regularization_a
    control%regularization_b = regularization_b
    control%max_lm_iterations = max_lm_iterations
    control%use_flux_objective = use_flux_objective
    control%gradient_tolerance = gradient_tolerance
    control%step_tolerance = step_tolerance
    control%lm_damping = lm_damping
    control%lm_damping_increase = lm_damping_increase
    control%lm_damping_decrease = lm_damping_decrease
    control%outer_boundary_density = outer_boundary_density

    print *, 'GORILLA_APPLETS: Loaded input data from global_transport_fit.inp'

end subroutine load_global_transport_fit_inp

end module global_transport_fit_settings_mod
