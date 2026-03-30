module km_benchmark_diagnostics_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    type collision_trace_t
        real(dp) :: tracer_mass = 0.0_dp
        real(dp) :: tracer_charge = 0.0_dp
        real(dp) :: energy_eV = 0.0_dp
        real(dp) :: total_collision_frequency_hz = 0.0_dp
        real(dp) :: total_efcolf = 0.0_dp
        real(dp) :: v0 = 0.0_dp
        real(dp), allocatable :: collision_frequency_hz(:)
        real(dp), allocatable :: efcolf(:)
        real(dp), allocatable :: enrat(:)
        real(dp), allocatable :: lambda(:)
        real(dp), allocatable :: velrat(:)
    end type collision_trace_t

    public :: collision_trace_t
    public :: build_log_energy_grid
    public :: compute_collision_trace
    public :: compute_energy_convolution
    public :: compute_maxwellian_energy_weights
    public :: compute_relative_change
    public :: compute_tracing_time_seconds
    public :: estimate_collision_time_seconds

contains

subroutine compute_collision_trace(tracer_mass, tracer_charge, masses, charges, &
    densities, temperatures, energy_eV, trace)

    use collis_ions, only: collis_init, lambda_alpha_beta
    use constants, only: echarge

    real(dp), intent(in) :: tracer_mass
    real(dp), intent(in) :: tracer_charge
    real(dp), intent(in) :: masses(:)
    real(dp), intent(in) :: charges(:)
    real(dp), intent(in) :: densities(:)
    real(dp), intent(in) :: temperatures(:)
    real(dp), intent(in) :: energy_eV
    type(collision_trace_t), intent(out) :: trace

    integer :: i, n_species

    n_species = size(masses)
    if (size(charges) /= n_species .or. size(densities) /= n_species .or. &
        size(temperatures) /= n_species) then
        error stop 'compute_collision_trace: inconsistent array sizes'
    end if

    if (allocated(trace%lambda)) deallocate(trace%lambda)
    if (allocated(trace%efcolf)) deallocate(trace%efcolf)
    if (allocated(trace%velrat)) deallocate(trace%velrat)
    if (allocated(trace%enrat)) deallocate(trace%enrat)
    if (allocated(trace%collision_frequency_hz)) deallocate(trace%collision_frequency_hz)
    allocate(trace%lambda(n_species))
    allocate(trace%efcolf(n_species))
    allocate(trace%velrat(n_species))
    allocate(trace%enrat(n_species))
    allocate(trace%collision_frequency_hz(n_species))

    call collis_init(tracer_mass, tracer_charge / echarge, masses, charges, &
        densities, temperatures, energy_eV, trace%v0, trace%efcolf, &
        trace%velrat, trace%enrat)

    do i = 1, n_species
        call lambda_alpha_beta(tracer_charge / echarge, charges(i), tracer_mass, &
            masses(i), energy_eV, temperatures(i), densities(i), densities(i), &
            trace%lambda(i))
    end do

    trace%tracer_mass = tracer_mass
    trace%tracer_charge = tracer_charge
    trace%energy_eV = energy_eV
    trace%collision_frequency_hz = trace%efcolf * trace%v0
    trace%total_efcolf = sum(trace%efcolf)
    trace%total_collision_frequency_hz = sum(trace%collision_frequency_hz)

end subroutine compute_collision_trace

real(dp) function estimate_collision_time_seconds(trace) result(collision_time_s)

    type(collision_trace_t), intent(in) :: trace

    collision_time_s = 0.0_dp
    if (trace%total_collision_frequency_hz > tiny(1.0_dp)) then
        collision_time_s = 1.0_dp / trace%total_collision_frequency_hz
    end if

end function estimate_collision_time_seconds

real(dp) function compute_tracing_time_seconds(override_time_s, trace_multiplier, &
    collision_time_s) result(trace_time_s)

    real(dp), intent(in) :: override_time_s
    real(dp), intent(in) :: trace_multiplier
    real(dp), intent(in) :: collision_time_s

    if (override_time_s > 0.0_dp) then
        trace_time_s = override_time_s
    else
        trace_time_s = trace_multiplier * collision_time_s
    end if

end function compute_tracing_time_seconds

subroutine build_log_energy_grid(center_energy_eV, min_factor, max_factor, &
    n_points, energy_grid_eV)

    real(dp), intent(in) :: center_energy_eV
    real(dp), intent(in) :: min_factor
    real(dp), intent(in) :: max_factor
    integer, intent(in) :: n_points
    real(dp), allocatable, intent(out) :: energy_grid_eV(:)

    integer :: i
    real(dp) :: log_max, log_min

    if (n_points < 1) error stop 'build_log_energy_grid: n_points must be positive'
    if (center_energy_eV <= 0.0_dp) error stop 'build_log_energy_grid: energy must be positive'
    if (min_factor <= 0.0_dp .or. max_factor <= 0.0_dp) then
        error stop 'build_log_energy_grid: factors must be positive'
    end if

    allocate(energy_grid_eV(n_points))
    if (n_points == 1) then
        energy_grid_eV(1) = center_energy_eV
        return
    end if

    log_min = log(center_energy_eV * min_factor)
    log_max = log(center_energy_eV * max_factor)
    do i = 1, n_points
        energy_grid_eV(i) = exp(log_min + real(i - 1, dp) * &
            (log_max - log_min) / real(n_points - 1, dp))
    end do

end subroutine build_log_energy_grid

subroutine compute_maxwellian_energy_weights(energy_eV, temperature_eV, weights)

    real(dp), intent(in) :: energy_eV(:)
    real(dp), intent(in) :: temperature_eV
    real(dp), allocatable, intent(out) :: weights(:)

    real(dp) :: normalization

    if (temperature_eV <= 0.0_dp) then
        error stop 'compute_maxwellian_energy_weights: temperature must be positive'
    end if

    allocate(weights(size(energy_eV)))
    ! libneo uses a generalized Gauss-Laguerre weight with alpha=5/2 for D11.
    weights = (max(energy_eV, 0.0_dp) / temperature_eV)**2.5_dp * &
        exp(-energy_eV / temperature_eV)
    normalization = sum(weights)
    if (normalization <= tiny(1.0_dp)) then
        error stop 'compute_maxwellian_energy_weights: invalid normalization'
    end if
    weights = weights / normalization

end subroutine compute_maxwellian_energy_weights

subroutine compute_energy_convolution(energy_eV, transport_values, &
    temperature_eV, convolved_value, weights)

    real(dp), intent(in) :: energy_eV(:)
    real(dp), intent(in) :: transport_values(:)
    real(dp), intent(in) :: temperature_eV
    real(dp), intent(out) :: convolved_value
    real(dp), allocatable, intent(out), optional :: weights(:)

    real(dp), allocatable :: energy_weights(:)
    real(dp) :: denominator, numerator
    integer :: i

    if (size(energy_eV) /= size(transport_values)) then
        error stop 'compute_energy_convolution: inconsistent array sizes'
    end if

    if (size(energy_eV) == 1) then
        convolved_value = transport_values(1)
        if (present(weights)) then
            allocate(weights(1))
            weights(1) = 1.0_dp
        end if
        return
    end if

    call compute_maxwellian_energy_weights(energy_eV, temperature_eV, energy_weights)

    numerator = 0.0_dp
    denominator = 0.0_dp
    do i = 1, size(energy_eV) - 1
        numerator = numerator + 0.5_dp * ( &
            energy_weights(i) * transport_values(i) + &
            energy_weights(i + 1) * transport_values(i + 1)) * &
            (energy_eV(i + 1) - energy_eV(i))
        denominator = denominator + 0.5_dp * ( &
            energy_weights(i) + energy_weights(i + 1)) * &
            (energy_eV(i + 1) - energy_eV(i))
    end do

    if (denominator <= tiny(1.0_dp)) then
        error stop 'compute_energy_convolution: invalid denominator'
    end if

    convolved_value = numerator / denominator
    if (present(weights)) call move_alloc(energy_weights, weights)

end subroutine compute_energy_convolution

real(dp) function compute_relative_change(previous_value, current_value) result(relative_change)

    real(dp), intent(in) :: previous_value
    real(dp), intent(in) :: current_value
    real(dp) :: scale

    scale = max(max(abs(previous_value), abs(current_value)), tiny(1.0_dp))
    relative_change = abs(current_value - previous_value) / scale

end function compute_relative_change

end module km_benchmark_diagnostics_mod
