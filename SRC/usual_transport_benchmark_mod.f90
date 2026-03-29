module usual_transport_benchmark_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    type, public :: usual_transport_result_t
        integer :: lost_particles = 0
        integer :: target_boundary_index = 0
        real(dp) :: density = 0.0_dp
        real(dp) :: density_log_gradient_per_s = 0.0_dp
        real(dp) :: density_profile_reference_s = 0.0_dp
        real(dp) :: surface_s = 0.0_dp
        real(dp) :: target_area = 0.0_dp
        real(dp) :: target_boundary_s = 0.0_dp
        real(dp) :: temperature_eV = 0.0_dp
        real(dp) :: total_time = 0.0_dp
        real(dp) :: v_E = 0.0_dp
        real(dp) :: weighted_crossings = 0.0_dp
        real(dp) :: weighted_flux_density = 0.0_dp
    end type usual_transport_result_t

    public :: calc_usual_transport_benchmark
    public :: run_usual_transport_benchmark_case

contains

subroutine calc_usual_transport_benchmark()

    use usual_transport_benchmark_settings_mod, only: boole_precalc_collisions, collision_operator, density, &
        density_log_gradient_per_s, density_profile_reference_s, electron_density, electron_temperature_eV, &
        filename_boundary_flux, filename_transport_summary, i_integrator_type, ion_density, ion_temperature_eV, &
        load_usual_transport_benchmark_inp, n_particles, surface_s, temperature_eV, total_time, tracer_species, v_E

    type(usual_transport_result_t) :: result

    call load_usual_transport_benchmark_inp()
    call run_usual_transport_benchmark_case(surface_s, tracer_species, density, electron_density, electron_temperature_eV, &
        ion_density, ion_temperature_eV, temperature_eV, total_time, n_particles, density_log_gradient_per_s, &
        density_profile_reference_s, v_E, collision_operator, boole_precalc_collisions, i_integrator_type, result, .true., &
        filename_transport_summary, filename_boundary_flux)

    print *, 'Usual transport benchmark complete.'
    print *, 'surface_s = ', result%surface_s
    print *, 'target_boundary = ', result%target_boundary_s
    print *, 'weighted_flux_density [1/cm^2/s] = ', result%weighted_flux_density
    print *, 'lost_particles = ', result%lost_particles

end subroutine calc_usual_transport_benchmark

subroutine run_usual_transport_benchmark_case(surface_s, tracer_species, density, electron_density, electron_temperature_eV, &
    ion_density, ion_temperature_eV, temperature_eV, total_time, n_particles, density_log_gradient_per_s, &
    density_profile_reference_s, v_E, collision_operator, boole_precalc_collisions, i_integrator_type, result, &
    write_outputs, summary_filename, boundary_filename)

    use field_mod, only: ipert
    use gorilla_settings_mod, only: coord_system
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use tetra_grid_mod, only: tetra_grid
    use tetra_physics_mod, only: make_tetra_physics
    use transport_benchmark_utils_mod, only: compute_boundary_areas, prepare_minimal_transport_output, &
        select_target_boundary, set_active_species_parameters, set_flux_surface_constant_eps_phi_from_vE
    use utils_data_pre_and_post_processing_mod, only: get_ipert, set_seed_for_random_numbers
    use utils_self_consistent_ef_mod, only: associate_flux_labels_with_tetrahedra_and_vertices, calc_s_shell_volumes, &
        print_errors_for_bad_inputs
    use volume_integrals_and_sqrt_g_mod, only: calc_volume_integrals_in_flux_coordinates

    real(dp), intent(in) :: density
    real(dp), intent(in) :: density_log_gradient_per_s
    real(dp), intent(in) :: density_profile_reference_s
    real(dp), intent(in) :: electron_density
    real(dp), intent(in) :: electron_temperature_eV
    real(dp), intent(in) :: ion_density
    real(dp), intent(in) :: ion_temperature_eV
    real(dp), intent(in) :: surface_s
    real(dp), intent(in) :: temperature_eV
    real(dp), intent(in) :: total_time
    real(dp), intent(in) :: v_E
    integer, intent(in) :: collision_operator
    integer, intent(in) :: i_integrator_type
    integer, intent(in) :: n_particles
    integer, intent(in) :: tracer_species
    logical, intent(in) :: boole_precalc_collisions
    type(usual_transport_result_t), intent(out) :: result
    logical, intent(in), optional :: write_outputs
    character(len=*), intent(in), optional :: boundary_filename
    character(len=*), intent(in), optional :: summary_filename

    real(dp), allocatable :: boundary_area(:)
    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: shell_volumes(:)
    real(dp), dimension(3) :: x_start
    logical :: reuse_initialized_grid
    integer :: target_boundary_index

    call set_seed_for_random_numbers()
    call configure_benchmark_input(surface_s, tracer_species, density, electron_density, electron_temperature_eV, &
        ion_density, ion_temperature_eV, temperature_eV, total_time, n_particles, density_log_gradient_per_s, &
        density_profile_reference_s, collision_operator, boole_precalc_collisions, i_integrator_type)
    reuse_initialized_grid = allocated(tetra_grid)
    call get_ipert()

    call set_active_species_parameters(tracer_species)
    if (.not.reuse_initialized_grid) call initialize_gorilla(1, ipert)
    x_start = (/surface_s, 0.0_dp, 0.0_dp/)
    call set_flux_surface_constant_eps_phi_from_vE(x_start, temperature_eV, v_E)

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

    call prepare_boundary_geometry(boundary_s, shell_volumes)
    call select_target_boundary(boundary_s, surface_s, target_boundary_index)
    call compute_boundary_areas(boundary_s, shell_volumes, boundary_area)
    call run_boundary_flux_transport_estimate(boundary_s, target_boundary_index, n_particles, result%lost_particles)
    call fill_transport_result(surface_s, temperature_eV, v_E, density, density_log_gradient_per_s, &
        density_profile_reference_s, total_time, n_particles, target_boundary_index, boundary_s, boundary_area, result)

    if (present(write_outputs)) then
        if (write_outputs) then
            if (.not.present(summary_filename)) then
                print *, 'Error: summary_filename missing for usual transport benchmark output.'
                stop
            end if
            if (.not.present(boundary_filename)) then
                print *, 'Error: boundary_filename missing for usual transport benchmark output.'
                stop
            end if
            call write_transport_summary(boundary_s, boundary_area, result, summary_filename, boundary_filename, n_particles)
        end if
    end if

end subroutine run_usual_transport_benchmark_case

subroutine configure_benchmark_input(surface_s, tracer_species, density, electron_density, electron_temperature_eV, &
    ion_density, ion_temperature_eV, temperature_eV, total_time, n_particles, density_log_gradient_per_s, &
    density_profile_reference_s, collision_operator, boole_precalc_collisions, i_integrator_type)

    use gorilla_applets_types_mod, only: in

    real(dp), intent(in) :: density
    real(dp), intent(in) :: density_log_gradient_per_s
    real(dp), intent(in) :: density_profile_reference_s
    real(dp), intent(in) :: electron_density
    real(dp), intent(in) :: electron_temperature_eV
    real(dp), intent(in) :: ion_density
    real(dp), intent(in) :: ion_temperature_eV
    real(dp), intent(in) :: surface_s
    real(dp), intent(in) :: temperature_eV
    real(dp), intent(in) :: total_time
    integer, intent(in) :: collision_operator
    integer, intent(in) :: i_integrator_type
    integer, intent(in) :: n_particles
    integer, intent(in) :: tracer_species
    logical, intent(in) :: boole_precalc_collisions

    in%boole_antithetic_variate = .true.
    in%boole_boltzmann_energies = .true.
    in%boole_calc_diffusion_coefficient = .false.
    in%boole_collisions = .true.
    in%boole_custom_background = .true.
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
    in%density = density
    in%density_log_gradient_per_s = density_log_gradient_per_s
    in%density_profile_reference_s = density_profile_reference_s
    in%energy_eV = temperature_eV
    in%background_density_cm3 = (/ion_density, electron_density/)
    in%background_temperature_eV = (/ion_temperature_eV, electron_temperature_eV/)
    in%i_integrator_type = i_integrator_type
    in%n_electric_potential_updates = 0
    in%n_particles = real(n_particles, dp)
    in%n_species = 1
    in%num_particles = n_particles
    in%time_step = total_time
    in%tracer_species = tracer_species
    in%update_dimension = 1

end subroutine configure_benchmark_input

subroutine prepare_boundary_geometry(boundary_s, shell_volumes)

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

end subroutine prepare_boundary_geometry

subroutine run_boundary_flux_transport_estimate(boundary_s, target_boundary_index, n_particles, lost_particles)

    use gorilla_applets_types_mod, only: exit_data, start
    use transport_benchmark_utils_mod, only: get_local_start_band
    use utils_data_pre_and_post_processing_mod, only: calc_collision_coefficients_for_all_tetrahedra, &
        initialize_exit_data, prepare_next_round_of_parallelised_particle_pushing
    use utils_self_consistent_ef_mod, only: allocate_start_type, parallelised_particle_pushing, &
        set_particle_type_specifications, set_rest_of_individual_particle_specifications, set_starting_positions

    real(dp), intent(in) :: boundary_s(:)
    integer, intent(in) :: n_particles
    integer, intent(in) :: target_boundary_index
    integer, intent(out) :: lost_particles

    real(dp), allocatable :: rand_matrix(:, :, :)
    real(dp) :: band_max_s
    real(dp) :: band_min_s

    allocate(rand_matrix(5, n_particles, 1))
    call random_number(rand_matrix)

    call allocate_start_type(n_particles)
    call set_particle_type_specifications()
    call initialize_exit_data(n_particles)
    call set_starting_positions(rand_matrix, species_in=(/1/))
    call get_local_start_band(boundary_s, target_boundary_index, band_min_s, band_max_s)
    start%x(1, :, 1) = band_min_s + (band_max_s - band_min_s) * rand_matrix(1, :, 1)
    call set_rest_of_individual_particle_specifications(rand_matrix, species_in=(/1/), n_particles_in=n_particles)
    call prepare_next_round_of_parallelised_particle_pushing(1)
    call calc_collision_coefficients_for_all_tetrahedra(1)
    call parallelised_particle_pushing(species=1, j=1, boole_diffusion_coefficient=.false., n_particles_in=n_particles)

    lost_particles = count(exit_data%lost(:, 1) /= 0)

end subroutine run_boundary_flux_transport_estimate

subroutine fill_transport_result(surface_s, temperature_eV, v_E, density, density_log_gradient_per_s, &
    density_profile_reference_s, total_time, n_particles, target_boundary_index, boundary_s, boundary_area, result)

    use gorilla_applets_types_mod, only: output

    real(dp), intent(in) :: boundary_area(:)
    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: density
    real(dp), intent(in) :: density_log_gradient_per_s
    real(dp), intent(in) :: density_profile_reference_s
    real(dp), intent(in) :: surface_s
    real(dp), intent(in) :: temperature_eV
    real(dp), intent(in) :: total_time
    real(dp), intent(in) :: v_E
    integer, intent(in) :: n_particles
    integer, intent(in) :: target_boundary_index
    type(usual_transport_result_t), intent(inout) :: result

    result%surface_s = surface_s
    result%target_boundary_index = target_boundary_index
    result%target_boundary_s = boundary_s(target_boundary_index)
    result%target_area = boundary_area(target_boundary_index)
    result%density = density
    result%temperature_eV = temperature_eV
    result%v_E = v_E
    result%density_log_gradient_per_s = density_log_gradient_per_s
    result%density_profile_reference_s = density_profile_reference_s
    result%total_time = total_time
    result%weighted_crossings = output%weighted_radial_flux(target_boundary_index)
    result%weighted_flux_density = result%weighted_crossings / (result%target_area * total_time * real(n_particles, dp))

end subroutine fill_transport_result

subroutine write_transport_summary(boundary_s, boundary_area, result, summary_filename, boundary_filename, n_particles)

    use gorilla_applets_types_mod, only: output

    real(dp), intent(in) :: boundary_area(:)
    real(dp), intent(in) :: boundary_s(:)
    integer, intent(in) :: n_particles
    type(usual_transport_result_t), intent(in) :: result
    character(len=*), intent(in) :: boundary_filename
    character(len=*), intent(in) :: summary_filename

    integer :: boundary_flux_unit
    integer :: i
    integer :: summary_unit
    real(dp) :: flux_density

    open(newunit=summary_unit, file=summary_filename, status='replace', action='write')
    write(summary_unit, '(A)') &
        'surface_s target_boundary_index target_boundary_s density_cm3 temperature_eV v_E density_log_gradient_per_s ' // &
        'density_profile_reference_s total_time_s target_area_cm2 weighted_flux_density_cm2s weighted_crossings lost_particles'
    write(summary_unit, '(12ES24.16,1X,I0)') result%surface_s, real(result%target_boundary_index, dp), result%target_boundary_s, &
        result%density, result%temperature_eV, result%v_E, result%density_log_gradient_per_s, &
        result%density_profile_reference_s, result%total_time, result%target_area, result%weighted_flux_density, &
        result%weighted_crossings, result%lost_particles
    close(summary_unit)

    open(newunit=boundary_flux_unit, file=boundary_filename, status='replace', action='write')
    write(boundary_flux_unit, '(A)') &
        'boundary_index boundary_s area_cm2 raw_crossings weighted_crossings weighted_flux_density_cm2s'
    do i = 1, size(boundary_s)
        flux_density = output%weighted_radial_flux(i) / (real(n_particles, dp) * result%total_time * boundary_area(i))
        write(boundary_flux_unit, '(I0,5(1X,ES24.16))') i, boundary_s(i), boundary_area(i), output%radial_flux(i), &
            output%weighted_radial_flux(i), flux_density
    end do
    close(boundary_flux_unit)

end subroutine write_transport_summary

end module usual_transport_benchmark_mod
