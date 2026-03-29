module global_transport_fit_types_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    public :: global_transport_experiment_t
    public :: global_transport_fit_control_t
    public :: global_transport_quality_t
    public :: global_transport_fit_result_t
    public :: global_transport_samples_t

    type global_transport_experiment_t
        real(dp), allocatable :: boundary_s(:)
        real(dp), allocatable :: shell_volumes(:)
        real(dp), allocatable :: source(:)
        real(dp), allocatable :: density(:)
        real(dp), allocatable :: density_variance(:)
        real(dp), allocatable :: flux(:)
        real(dp), allocatable :: flux_variance(:)
        real(dp), allocatable :: integrated_flux(:)
        real(dp), allocatable :: integrated_flux_variance(:)
    end type global_transport_experiment_t

    type global_transport_quality_t
        logical :: transport_supported = .false.
        integer :: n_supported_flux_boundaries = 0
        integer :: n_batches = 0
        integer :: n_trials = 0
        real(dp) :: density_source_relstd = 0.0d0
        real(dp) :: tracing_time = 0.0d0
        real(dp) :: lost_fraction = 0.0d0
    end type global_transport_quality_t

    type global_transport_fit_control_t
        integer :: n_knots_a = 5
        integer :: n_knots_b = 5
        integer :: max_lm_iterations = 30
        logical :: use_flux_objective = .true.
        real(dp) :: regularization_a = 1.0d-8
        real(dp) :: regularization_b = 1.0d-8
        real(dp) :: density_relative_sigma_floor = 1.0d-2
        real(dp) :: flux_relative_sigma_floor = 1.0d-2
        real(dp) :: density_support_fraction = 1.0d-3
        real(dp) :: lm_damping = 1.0d-4
        real(dp) :: lm_damping_increase = 10.0d0
        real(dp) :: lm_damping_decrease = 0.3d0
        real(dp) :: step_tolerance = 1.0d-10
        real(dp) :: gradient_tolerance = 1.0d-10
        real(dp) :: outer_boundary_density = 0.0d0
        real(dp) :: support_sigma_multiplier = 2.0d0
    end type global_transport_fit_control_t

    type global_transport_fit_result_t
        logical :: converged = .false.
        integer :: n_iterations = 0
        real(dp) :: objective = 0.0d0
        real(dp), allocatable :: parameters(:)
        real(dp), allocatable :: gradient(:)
        real(dp), allocatable :: covariance(:,:)
        real(dp), allocatable :: a_profile(:)
        real(dp), allocatable :: b_profile(:)
        real(dp), allocatable :: a_std(:)
        real(dp), allocatable :: b_std(:)
        real(dp), allocatable :: knot_s_a(:)
        real(dp), allocatable :: knot_s_b(:)
    end type global_transport_fit_result_t

    type global_transport_samples_t
        real(dp), allocatable :: source_weight_sum(:)
        real(dp), allocatable :: shell_time_sum(:)
        real(dp), allocatable :: boundary_weighted_flux_sum(:)
        real(dp) :: mean_tracing_time = 0.0d0
        integer :: n_particles = 0
        integer :: lost_particles = 0
    end type global_transport_samples_t

end module global_transport_fit_types_mod
