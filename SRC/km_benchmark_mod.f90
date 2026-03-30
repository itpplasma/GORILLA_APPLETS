module km_benchmark_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    public :: calc_km_benchmark

contains

subroutine calc_km_benchmark()

    use km_benchmark_settings_mod, only: background_charge, background_density, &
        background_mass_amu, background_temperature, collision_operator, &
        filename_output, i_integrator_type, load_km_benchmark_inp, &
        n_background_species, n_particles, n_surfaces, boole_precalc_collisions, &
        surface_s_values, temperature_eV, total_time, tracer_species, v_E
    use kramers_moyal_transport_mod, only: calc_km_d11_profile, km_d11_result_t

    type(km_d11_result_t) :: result
    integer, allocatable :: surface_indices(:)

    call load_km_benchmark_inp()
    call configure_km_input()
    call initialize_km_grid(tracer_species, temperature_eV, v_E)
    call map_s_values_to_indices(surface_s_values, n_surfaces, surface_indices)
    call calc_km_d11_profile(surface_indices, result)
    call write_km_csv(result, filename_output)

    print *, ''
    print *, 'KM benchmark complete. Output written to ', trim(filename_output)

end subroutine calc_km_benchmark

subroutine configure_km_input()

    use km_benchmark_settings_mod, only: background_charge, background_density, &
        background_mass_amu, background_temperature, collision_operator, &
        n_background_species, n_particles, boole_precalc_collisions, &
        temperature_eV, total_time, tracer_species, v_E, i_integrator_type
    use gorilla_applets_types_mod, only: in
    use utils_data_pre_and_post_processing_mod, only: set_custom_background
    use constants, only: amp

    integer :: n

    n = n_background_species

    in%boole_antithetic_variate = .false.
    in%boole_boltzmann_energies = .false.
    in%boole_calc_diffusion_coefficient = .true.
    in%boole_collisions = .true.
    in%boole_divertor_intersection = .false.
    in%boole_linear_density_simulation = .false.
    in%boole_linear_temperature_simulation = .false.
    in%boole_point_source = .false.
    in%boole_precalc_collisions = boole_precalc_collisions
    in%boole_preserve_energy_and_momentum_during_collisions = .false.
    in%boole_refined_sqrt_g = .false.
    in%boole_squared_moments = .false.
    in%boole_static_ne = .true.
    in%boole_write_boltzmann_density = .false.
    in%boole_write_electric_potential = .false.
    in%boole_write_exit_data = .false.
    in%boole_write_fourier_moments = .false.
    in%boole_write_grid_data = .false.
    in%boole_write_moments = .false.
    in%boole_write_prism_volumes = .false.
    in%boole_write_refined_prism_volumes = .false.
    in%boole_write_vertex_coordinates = .false.
    in%boole_write_vertex_indices = .false.
    in%collision_operator = collision_operator
    in%density = background_density(n)
    in%energy_eV = temperature_eV
    call set_custom_background(n, &
        background_density(1:n), &
        background_temperature(1:n), &
        background_mass_amu(1:n) * amp, &
        background_charge(1:n))
    in%i_integrator_type = i_integrator_type
    in%n_electric_potential_updates = 0
    in%n_particles = real(n_particles, dp)
    in%n_species = 1
    in%num_particles = n_particles
    in%time_step = total_time
    in%tracer_species = tracer_species
    in%update_dimension = 1

end subroutine configure_km_input

subroutine initialize_km_grid(tracer_species, temperature_eV, v_E)

    use field_mod, only: ipert
    use gorilla_settings_mod, only: coord_system
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use tetra_grid_mod, only: tetra_grid
    use tetra_physics_mod, only: make_tetra_physics
    use transport_benchmark_utils_mod, only: prepare_minimal_transport_output, &
        set_active_species_parameters, set_flux_surface_constant_eps_phi_from_vE
    use utils_data_pre_and_post_processing_mod, only: get_ipert, set_seed_for_random_numbers
    use utils_self_consistent_ef_mod, only: associate_flux_labels_with_tetrahedra_and_vertices, &
        calc_s_shell_volumes, print_errors_for_bad_inputs
    use volume_integrals_and_sqrt_g_mod, only: calc_volume_integrals_in_flux_coordinates

    integer, intent(in) :: tracer_species
    real(dp), intent(in) :: temperature_eV
    real(dp), intent(in) :: v_E

    real(dp), dimension(3) :: x_start
    logical :: reuse_initialized_grid

    call set_seed_for_random_numbers()
    reuse_initialized_grid = allocated(tetra_grid)
    call get_ipert()

    call set_active_species_parameters(tracer_species)
    if (.not. reuse_initialized_grid) call initialize_gorilla(1, ipert)
    if (abs(v_E) > 0.0_dp) then
        x_start = (/0.5_dp, 0.0_dp, 0.0_dp/)
        call set_flux_surface_constant_eps_phi_from_vE(x_start, temperature_eV, v_E)
    end if

    if (reuse_initialized_grid) then
        call make_tetra_physics(coord_system, ipert)
    else
        call initialize_gorilla(2, ipert)
        call set_active_species_parameters(tracer_species)
    end if
    call prepare_minimal_transport_output()
    call calc_volume_integrals_in_flux_coordinates()
    call associate_flux_labels_with_tetrahedra_and_vertices()
    call calc_s_shell_volumes()
    call print_errors_for_bad_inputs()

end subroutine initialize_km_grid

subroutine map_s_values_to_indices(surface_s_values, n_surfaces, surface_indices)

    use tetra_grid_mod, only: verts_sthetaphi
    use tetra_grid_settings_mod, only: grid_size

    real(dp), intent(in) :: surface_s_values(:)
    integer, intent(in) :: n_surfaces
    integer, allocatable, intent(out) :: surface_indices(:)

    real(dp), allocatable :: boundary_s(:)
    real(dp) :: best_dist, dist
    integer :: i, j, best_idx, n_valid, actual_n

    allocate(boundary_s(grid_size(1) + 1))
    do i = 1, grid_size(1) + 1
        boundary_s(i) = verts_sthetaphi(1, grid_size(3) * (i - 1) + 1)
    end do

    if (n_surfaces > 0) then
        actual_n = min(n_surfaces, size(surface_s_values))
    else
        actual_n = count(surface_s_values > 0.0_dp)
    end if

    if (actual_n == 0) then
        n_valid = grid_size(1) - 2
        allocate(surface_indices(n_valid))
        do i = 1, n_valid
            surface_indices(i) = i + 1
        end do
        print '(A,I0,A)', '  Using all ', n_valid, ' interior grid surfaces'
        return
    end if

    allocate(surface_indices(actual_n))
    do i = 1, actual_n
        best_idx = 2
        best_dist = abs(boundary_s(2) - surface_s_values(i))
        do j = 3, grid_size(1)
            dist = abs(boundary_s(j) - surface_s_values(i))
            if (dist < best_dist) then
                best_dist = dist
                best_idx = j
            end if
        end do
        surface_indices(i) = best_idx
        print '(A,F8.5,A,I0,A,F8.5)', '  Mapped s=', surface_s_values(i), &
            ' -> index ', best_idx, ' (s_grid=', boundary_s(best_idx), ')'
    end do

end subroutine map_s_values_to_indices

subroutine write_km_csv(result, filename)

    use kramers_moyal_transport_mod, only: km_d11_result_t

    type(km_d11_result_t), intent(in) :: result
    character(len=*), intent(in) :: filename

    integer :: csv_unit, i

    open(newunit=csv_unit, file=filename, status='replace', action='write')
    write(csv_unit, '(A)') 's,d11_s,convection_A'
    do i = 1, result%n_surfaces
        write(csv_unit, '(ES24.16,A,ES24.16,A,ES24.16)') &
            result%boundary_s(result%surface_indices(i)), ',', &
            result%d11(i), ',', result%convection(i)
    end do
    close(csv_unit)

end subroutine write_km_csv

end module km_benchmark_mod
