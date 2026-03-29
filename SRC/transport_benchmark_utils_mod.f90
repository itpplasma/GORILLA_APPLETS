module transport_benchmark_utils_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use global_transport_fit_math_mod, only: compute_geometry_from_boundaries

    implicit none

    private

    public :: fit_transport_coefficients
    public :: set_active_species_parameters
    public :: set_flux_surface_constant_eps_phi_from_vE
    public :: get_species_parameters
    public :: prepare_minimal_transport_output
    public :: select_target_boundary
    public :: compute_boundary_areas
    public :: get_local_start_band
    public :: recover_transport_coefficients_from_flux_pair

contains

subroutine fit_transport_coefficients(time, delta_s, delta_s_squared, A, B, fit_start_fraction, fit_end_fraction)

    use llsq_mod, only: llsq

    real(dp), dimension(:), intent(in) :: time
    real(dp), dimension(:), intent(in) :: delta_s
    real(dp), dimension(:), intent(in) :: delta_s_squared
    real(dp), intent(out) :: A
    real(dp), intent(out) :: B
    real(dp), intent(in), optional :: fit_end_fraction
    real(dp), intent(in), optional :: fit_start_fraction

    real(dp), dimension(:), allocatable :: diffusion_data
    real(dp) :: effective_fit_end_fraction
    real(dp) :: effective_fit_start_fraction
    real(dp) :: offset
    integer :: first_fit_index
    integer :: last_fit_index
    integer :: n_samples

    n_samples = size(time)
    effective_fit_start_fraction = 0.2_dp
    effective_fit_end_fraction = 1.0_dp
    if (present(fit_start_fraction)) effective_fit_start_fraction = fit_start_fraction
    if (present(fit_end_fraction)) effective_fit_end_fraction = fit_end_fraction
    effective_fit_start_fraction = max(0.0_dp, min(effective_fit_start_fraction, 0.999_dp))
    effective_fit_end_fraction = max(effective_fit_start_fraction + 1.0e-3_dp, min(effective_fit_end_fraction, 1.0_dp))
    first_fit_index = max(1, ceiling(real(n_samples, dp) * effective_fit_start_fraction))
    last_fit_index = min(n_samples, ceiling(real(n_samples, dp) * effective_fit_end_fraction))
    if (last_fit_index <= first_fit_index) last_fit_index = min(n_samples, first_fit_index + 1)
    allocate(diffusion_data(n_samples))

    call llsq(int(last_fit_index - first_fit_index + 1, kind=8), time(first_fit_index:last_fit_index), &
              delta_s(first_fit_index:last_fit_index), A, offset)

    diffusion_data = delta_s_squared - 2.0_dp * A * time * delta_s + A**2 * time**2

    call llsq(int(last_fit_index - first_fit_index + 1, kind=8), time(first_fit_index:last_fit_index), &
              diffusion_data(first_fit_index:last_fit_index), B, offset)
    B = 0.5_dp * B

    deallocate(diffusion_data)

end subroutine fit_transport_coefficients

subroutine set_active_species_parameters(species_id)

    use tetra_physics_mod, only: cm_over_e, particle_charge, particle_mass

    integer, intent(in) :: species_id
    real(dp) :: charge
    real(dp) :: cm_over_e_local
    real(dp) :: mass

    call get_species_parameters(species_id, charge, mass, cm_over_e_local)

    particle_charge = charge
    particle_mass = mass
    cm_over_e = cm_over_e_local

end subroutine set_active_species_parameters

subroutine get_species_parameters(species_id, charge, mass, cm_over_e)

    use constants, only: ame, amp, clight, echarge

    integer, intent(in) :: species_id
    real(dp), intent(out) :: charge
    real(dp), intent(out) :: cm_over_e
    real(dp), intent(out) :: mass

    select case (species_id)
        case (1)
            charge = -echarge
            mass = ame
            cm_over_e = -clight * ame / echarge
        case (2)
            charge = echarge
            mass = 2.0_dp * amp
            cm_over_e = 2.0_dp * clight * amp / echarge
        case (3)
            charge = 2.0_dp * echarge
            mass = 4.0_dp * amp
            cm_over_e = 2.0_dp * clight * amp / echarge
        case (4)
            charge = 74.0_dp * echarge
            mass = 184.0_dp * amp
            cm_over_e = 184.0_dp * clight * amp / (74.0_dp * echarge)
        case default
            print *, 'Error: unsupported tracer_species ', species_id
            stop
    end select

end subroutine get_species_parameters

subroutine set_flux_surface_constant_eps_phi_from_vE(x_start, reference_energy_eV, normalized_vE)

    use constants, only: clight, ev2erg, pi
    use field_divB0_mod, only: field
    use gorilla_settings_mod, only: set_eps_Phi
    use magdata_in_symfluxcoor_mod, only: psitor_max
    use splint_vmec_data_mod, only: splint_vmec_data
    use sub_alpha_lifetime_can_mod, only: integrate_mfl_can
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_physics_mod, only: mag_axis_R0, mag_axis_Z0, particle_mass
    use vector_potential_mod, only: torflux

    real(dp), dimension(3), intent(in) :: x_start
    real(dp), intent(in) :: normalized_vE
    real(dp), intent(in) :: reference_energy_eV

    real(dp) :: aiota
    real(dp) :: alam
    real(dp) :: A_phi
    real(dp) :: A_theta
    real(dp) :: bmod00
    real(dp) :: B_p
    real(dp) :: B_r
    real(dp) :: B_z
    real(dp) :: dBrdR
    real(dp) :: dBrdZ
    real(dp) :: dBrdp
    real(dp) :: dBpdR
    real(dp) :: dBpdZ
    real(dp) :: dBpdp
    real(dp) :: dBzdR
    real(dp) :: dBzdZ
    real(dp) :: dBzdp
    real(dp) :: dA_phi_ds
    real(dp) :: dA_theta_ds
    real(dp) :: dl_dp
    real(dp) :: dl_ds
    real(dp) :: dl_dt
    real(dp) :: dphi
    real(dp) :: dR_dp
    real(dp) :: dR_ds
    real(dp) :: dR_dt
    real(dp) :: ds_dr
    real(dp) :: dZ_dp
    real(dp) :: dZ_ds
    real(dp) :: dZ_dt
    real(dp) :: eps_Phi
    integer :: ierr
    integer, parameter :: npoi = 100000
    real(dp) :: phibeg
    real(dp) :: rbeg
    real(dp) :: s
    real(dp) :: theta
    real(dp) :: varphi
    real(dp) :: vmod
    real(dp), dimension(:), allocatable :: bstart
    real(dp), dimension(:), allocatable :: volstart
    real(dp), dimension(:, :), allocatable :: xstart
    real(dp) :: zbeg
    real(dp) :: Z
    real(dp) :: R

    vmod = sqrt(2.0_dp * reference_energy_eV * ev2erg / particle_mass)

    s = x_start(1)
    select case (grid_kind)
        case (2)
            call field(mag_axis_R0, 0.0_dp, mag_axis_Z0, B_r, B_p, B_z, dBrdR, dBrdp, dBrdZ, &
                       dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
            bmod00 = sqrt(B_r**2 + B_p**2 + B_z**2)
            dA_theta_ds = psitor_max
            ds_dr = 2.0_dp * sqrt(s) * sqrt(pi * bmod00 / psitor_max)
        case (3)
            allocate(xstart(3, npoi))
            allocate(bstart(npoi))
            allocate(volstart(npoi))

            dphi = 2.0_dp * pi / 100.0_dp
            rbeg = x_start(1)
            phibeg = x_start(3)
            zbeg = x_start(2)

            call integrate_mfl_can(npoi, dphi, rbeg, phibeg, zbeg, xstart, bstart, volstart, bmod00, ierr)

            theta = x_start(2)
            varphi = x_start(3)
            call splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, &
                                  dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
            ds_dr = 2.0_dp * sqrt(s) * sqrt(pi * bmod00 / torflux)

            deallocate(xstart)
            deallocate(bstart)
            deallocate(volstart)
        case default
            print *, 'Error: normalized_vE benchmark is only implemented for flux grids.'
            stop
    end select

    eps_Phi = normalized_vE * vmod * bmod00 / (dA_theta_ds * ds_dr * clight)
    call set_eps_Phi(eps_Phi)

end subroutine set_flux_surface_constant_eps_phi_from_vE

subroutine prepare_minimal_transport_output()

    use gorilla_applets_types_mod, only: in, moment_specs, output
    use tetra_grid_mod, only: ntetr
    use tetra_grid_settings_mod, only: grid_size
    use utils_data_pre_and_post_processing_mod, only: deallocate_output

    integer :: n_prisms

    n_prisms = ntetr / 3

    moment_specs%boole_squared_moments = .false.
    moment_specs%moments_selector = 0
    moment_specs%n_moments = 0
    moment_specs%n_triangles = 0
    moment_specs%n_fourier_modes = 0

    call deallocate_output()

    allocate(output%prism_volumes(n_prisms))
    output%prism_volumes = 0.0_dp

    allocate(output%radial_flux(grid_size(1) + 1))
    output%radial_flux = 0.0_dp
    allocate(output%weighted_radial_flux(grid_size(1) + 1))
    output%weighted_radial_flux = 0.0_dp

    if (in%boole_refined_sqrt_g) then
        allocate(output%refined_prism_volumes(n_prisms))
        output%refined_prism_volumes = 0.0_dp
    end if

end subroutine prepare_minimal_transport_output

subroutine select_target_boundary(boundary_s, target_s, target_index)

    real(dp), dimension(:), intent(in) :: boundary_s
    real(dp), intent(in) :: target_s
    integer, intent(out) :: target_index

    integer :: i
    real(dp) :: best_distance
    real(dp) :: distance

    target_index = 1
    best_distance = abs(boundary_s(1) - target_s)

    do i = 2, size(boundary_s)
        distance = abs(boundary_s(i) - target_s)
        if (distance < best_distance) then
            best_distance = distance
            target_index = i
        end if
    end do

end subroutine select_target_boundary

subroutine compute_boundary_areas(boundary_s, shell_volumes, boundary_area)

    real(dp), dimension(:), intent(in) :: boundary_s
    real(dp), dimension(:), intent(in) :: shell_volumes
    real(dp), dimension(:), allocatable, intent(out) :: boundary_area
    real(dp), allocatable :: cell_centers(:)

    call compute_geometry_from_boundaries(boundary_s, shell_volumes, cell_centers, boundary_area)

end subroutine compute_boundary_areas

subroutine get_local_start_band(boundary_s, target_index, band_min_s, band_max_s)

    real(dp), dimension(:), intent(in) :: boundary_s
    integer, intent(in) :: target_index
    real(dp), intent(out) :: band_min_s
    real(dp), intent(out) :: band_max_s

    band_min_s = boundary_s(max(1, target_index - 1))
    band_max_s = boundary_s(min(size(boundary_s), target_index + 1))

end subroutine get_local_start_band

subroutine recover_transport_coefficients_from_flux_pair(flux_density_1, flux_density_2, density_1, density_2, &
    density_gradient_1, density_gradient_2, a_coeff, b_coeff)

    real(dp), intent(in) :: flux_density_1
    real(dp), intent(in) :: flux_density_2
    real(dp), intent(in) :: density_1
    real(dp), intent(in) :: density_2
    real(dp), intent(in) :: density_gradient_1
    real(dp), intent(in) :: density_gradient_2
    real(dp), intent(out) :: a_coeff
    real(dp), intent(out) :: b_coeff

    real(dp) :: normalized_flux_1
    real(dp) :: normalized_flux_2
    real(dp) :: gradient_delta

    gradient_delta = density_gradient_2 - density_gradient_1
    if (abs(gradient_delta) < 1.0d-12) then
        print *, 'Error: density gradients must differ to recover local transport coefficients.'
        stop
    end if

    normalized_flux_1 = flux_density_1 / max(density_1, 1.0d-30)
    normalized_flux_2 = flux_density_2 / max(density_2, 1.0d-30)

    b_coeff = (normalized_flux_1 - normalized_flux_2) / gradient_delta
    a_coeff = normalized_flux_1 + b_coeff * density_gradient_1

end subroutine recover_transport_coefficients_from_flux_pair

end module transport_benchmark_utils_mod
