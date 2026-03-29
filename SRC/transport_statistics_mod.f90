module transport_statistics_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan

    implicit none

    private

    public :: build_gradient_grid
    public :: compute_density_source_relstd
    public :: compute_mean_and_variance
    public :: build_density_support_mask
    public :: build_inverse_variance_weights
    public :: build_supported_boundary_mask
    public :: count_supported_values
    public :: count_supported_boundaries
    public :: fit_linear_transport_response
    public :: sanitize_transport_output_triplet
    public :: transport_signal_supported

contains

subroutine build_gradient_grid(gradient_min, gradient_max, n_gradients, gradients)

    real(dp), intent(in) :: gradient_max
    real(dp), intent(in) :: gradient_min
    integer, intent(in) :: n_gradients
    real(dp), allocatable, intent(out) :: gradients(:)

    integer :: i

    allocate(gradients(n_gradients))
    if (n_gradients == 1) then
        gradients(1) = 0.5_dp * (gradient_min + gradient_max)
        return
    end if

    do i = 1, n_gradients
        gradients(i) = gradient_min + real(i - 1, dp) * (gradient_max - gradient_min) / real(n_gradients - 1, dp)
    end do

end subroutine build_gradient_grid

subroutine compute_mean_and_variance(samples, mean_values, variance_of_mean)

    real(dp), intent(in) :: samples(:, :)
    real(dp), allocatable, intent(out) :: mean_values(:)
    real(dp), allocatable, intent(out) :: variance_of_mean(:)

    integer :: i_batch
    integer :: n_batches
    real(dp), allocatable :: centered(:, :)

    n_batches = size(samples, 2)
    allocate(mean_values(size(samples, 1)))
    allocate(variance_of_mean(size(samples, 1)))

    mean_values = sum(samples, dim=2) / real(max(n_batches, 1), dp)
    if (n_batches <= 1) then
        variance_of_mean = 1.0d-30
        return
    end if

    allocate(centered(size(samples, 1), n_batches))
    do i_batch = 1, n_batches
        centered(:, i_batch) = samples(:, i_batch) - mean_values
    end do
    variance_of_mean = sum(centered**2, dim=2) / real(n_batches * (n_batches - 1), dp)
    variance_of_mean = max(variance_of_mean, 1.0d-30)

end subroutine compute_mean_and_variance

subroutine compute_density_source_relstd(source_density, shell_density, relstd)

    real(dp), intent(in) :: shell_density(:)
    real(dp), intent(in) :: source_density(:)
    real(dp), intent(out) :: relstd

    integer :: i
    integer :: n_selected
    real(dp) :: mean_ratio
    real(dp) :: max_source
    real(dp) :: variance_ratio
    real(dp), allocatable :: ratio(:)

    max_source = maxval(source_density)
    if (max_source <= 0.0_dp) then
        relstd = 0.0_dp
        return
    end if

    allocate(ratio(size(source_density)))
    n_selected = 0
    do i = 1, size(source_density)
        if (source_density(i) > 1.0d-6 * max_source) then
            n_selected = n_selected + 1
            ratio(n_selected) = shell_density(i) / source_density(i)
        end if
    end do
    if (n_selected <= 1) then
        relstd = 0.0_dp
        return
    end if

    mean_ratio = sum(ratio(1:n_selected)) / real(n_selected, dp)
    if (abs(mean_ratio) <= 1.0d-30) then
        relstd = 0.0_dp
        return
    end if
    variance_ratio = sum((ratio(1:n_selected) - mean_ratio)**2) / real(n_selected - 1, dp)
    relstd = sqrt(max(variance_ratio, 0.0_dp)) / abs(mean_ratio)

end subroutine compute_density_source_relstd

subroutine build_density_support_mask(source_density, shell_density, support_fraction, supported_mask)

    real(dp), intent(in) :: shell_density(:)
    real(dp), intent(in) :: source_density(:)
    real(dp), intent(in) :: support_fraction
    logical, allocatable, intent(out) :: supported_mask(:)

    real(dp) :: density_scale
    real(dp) :: density_threshold
    real(dp) :: source_scale
    real(dp) :: source_threshold

    allocate(supported_mask(size(shell_density)))
    source_scale = max(maxval(abs(source_density)), 1.0d-30)
    density_scale = max(maxval(abs(shell_density)), 1.0d-30)
    source_threshold = max(support_fraction, 0.0_dp) * source_scale
    density_threshold = max(support_fraction, 0.0_dp) * density_scale
    supported_mask = abs(source_density) >= source_threshold .or. abs(shell_density) >= density_threshold

end subroutine build_density_support_mask

subroutine build_inverse_variance_weights(values, variance, sigma_floor_fraction, supported_mask, weights)

    real(dp), intent(in) :: values(:)
    real(dp), intent(in) :: variance(:)
    real(dp), intent(in) :: sigma_floor_fraction
    logical, intent(in) :: supported_mask(:)
    real(dp), allocatable, intent(out) :: weights(:)

    integer :: i
    real(dp) :: sigma_floor
    real(dp) :: value_scale

    allocate(weights(size(values)))
    weights = 0.0_dp
    if (any(supported_mask)) then
        value_scale = max(maxval(abs(values), mask=supported_mask), 1.0d-30)
    else
        value_scale = 1.0d0
    end if
    sigma_floor = max(max(sigma_floor_fraction, 0.0_dp) * value_scale, 1.0d-15)

    do i = 1, size(values)
        if (.not. supported_mask(i)) cycle
        weights(i) = 1.0_dp / max(variance(i), sigma_floor**2)
    end do

end subroutine build_inverse_variance_weights

subroutine build_supported_boundary_mask(boundary_flux, boundary_flux_variance, sigma_multiplier, supported_mask)

    real(dp), intent(in) :: boundary_flux(:)
    real(dp), intent(in) :: boundary_flux_variance(:)
    real(dp), intent(in) :: sigma_multiplier
    logical, allocatable, intent(out) :: supported_mask(:)

    integer :: i

    allocate(supported_mask(size(boundary_flux)))
    supported_mask = .false.

    do i = 2, size(boundary_flux)
        supported_mask(i) = abs(boundary_flux(i)) > sigma_multiplier * sqrt(max(boundary_flux_variance(i), 0.0_dp))
    end do

end subroutine build_supported_boundary_mask

integer function count_supported_values(values, value_variance, sigma_multiplier) result(n_supported)

    real(dp), intent(in) :: values(:)
    real(dp), intent(in) :: value_variance(:)
    real(dp), intent(in) :: sigma_multiplier

    n_supported = count(abs(values) > sigma_multiplier * sqrt(max(value_variance, 0.0_dp)))

end function count_supported_values

integer function count_supported_boundaries(boundary_flux, boundary_flux_variance, sigma_multiplier) result(n_supported)

    real(dp), intent(in) :: boundary_flux(:)
    real(dp), intent(in) :: boundary_flux_variance(:)
    real(dp), intent(in) :: sigma_multiplier

    logical, allocatable :: supported_mask(:)

    call build_supported_boundary_mask(boundary_flux, boundary_flux_variance, sigma_multiplier, supported_mask)
    n_supported = count(supported_mask)

end function count_supported_boundaries

subroutine fit_linear_transport_response(gradients, normalized_flux, normalized_flux_variance, a_coeff, b_coeff, a_std, b_std, &
    valid_fit)

    real(dp), intent(in) :: gradients(:)
    real(dp), intent(in) :: normalized_flux(:)
    real(dp), intent(in) :: normalized_flux_variance(:)
    real(dp), intent(out) :: a_coeff
    real(dp), intent(out) :: b_coeff
    real(dp), intent(out) :: a_std
    real(dp), intent(out) :: b_std
    logical, intent(out) :: valid_fit

    real(dp) :: determinant
    real(dp) :: rhs_1
    real(dp) :: rhs_2
    real(dp) :: s00
    real(dp) :: s01
    real(dp) :: s11
    real(dp) :: weight
    real(dp) :: x_value
    integer :: i

    s00 = 0.0_dp
    s01 = 0.0_dp
    s11 = 0.0_dp
    rhs_1 = 0.0_dp
    rhs_2 = 0.0_dp

    do i = 1, size(gradients)
        if (normalized_flux_variance(i) <= 0.0_dp) cycle
        weight = 1.0_dp / normalized_flux_variance(i)
        x_value = -gradients(i)
        s00 = s00 + weight
        s01 = s01 + weight * x_value
        s11 = s11 + weight * x_value * x_value
        rhs_1 = rhs_1 + weight * normalized_flux(i)
        rhs_2 = rhs_2 + weight * x_value * normalized_flux(i)
    end do

    determinant = s00 * s11 - s01 * s01
    if (determinant <= 1.0d-20) then
        a_coeff = 0.0_dp
        b_coeff = 0.0_dp
        a_std = huge(1.0_dp)
        b_std = huge(1.0_dp)
        valid_fit = .false.
        return
    end if

    a_coeff = (rhs_1 * s11 - rhs_2 * s01) / determinant
    b_coeff = (s00 * rhs_2 - s01 * rhs_1) / determinant
    a_std = sqrt(max(s11 / determinant, 0.0_dp))
    b_std = sqrt(max(s00 / determinant, 0.0_dp))
    valid_fit = .true.

end subroutine fit_linear_transport_response

subroutine sanitize_transport_output_triplet(valid_value, mean_value, std_value, sanitized_mean, sanitized_std, sanitized_2sigma)

    logical, intent(in) :: valid_value
    real(dp), intent(in) :: mean_value
    real(dp), intent(in) :: std_value
    real(dp), intent(out) :: sanitized_mean
    real(dp), intent(out) :: sanitized_std
    real(dp), intent(out) :: sanitized_2sigma

    real(dp) :: max_std
    real(dp) :: quiet_nan_value

    quiet_nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    if (.not. valid_value) then
        sanitized_mean = quiet_nan_value
        sanitized_std = quiet_nan_value
        sanitized_2sigma = quiet_nan_value
        return
    end if

    if (.not. ieee_is_finite(mean_value)) then
        sanitized_mean = quiet_nan_value
    else
        sanitized_mean = mean_value
    end if

    if (.not. ieee_is_finite(std_value) .or. std_value < 0.0_dp) then
        sanitized_std = quiet_nan_value
        sanitized_2sigma = quiet_nan_value
        return
    end if

    max_std = 0.5_dp * huge(1.0_dp)
    sanitized_std = min(std_value, max_std)
    sanitized_2sigma = 2.0_dp * sanitized_std

end subroutine sanitize_transport_output_triplet

logical function transport_signal_supported(n_supported_flux_boundaries, min_supported_flux_boundaries, density_source_relstd, &
    min_density_source_relstd) result(is_supported)

    integer, intent(in) :: min_supported_flux_boundaries
    integer, intent(in) :: n_supported_flux_boundaries
    real(dp), intent(in) :: density_source_relstd
    real(dp), intent(in) :: min_density_source_relstd

    is_supported = n_supported_flux_boundaries >= min_supported_flux_boundaries
    is_supported = is_supported .and. density_source_relstd >= min_density_source_relstd

end function transport_signal_supported

end module transport_statistics_mod
