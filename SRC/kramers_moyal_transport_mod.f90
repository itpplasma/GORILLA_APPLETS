module kramers_moyal_transport_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    public :: calc_km_d11_profile
    public :: km_d11_result_t
    public :: set_trace_time_for_species
    public :: write_km_csv

    type km_d11_result_t
        real(dp), allocatable :: boundary_s(:)
        real(dp), allocatable :: collision_time_s(:)
        real(dp), allocatable :: convolved_d11(:)
        real(dp), allocatable :: convection(:)
        real(dp), allocatable :: d11(:)
        real(dp), allocatable :: energy_eV(:)
        real(dp), allocatable :: total_collision_frequency_hz(:)
        real(dp), allocatable :: trace_time_s(:)
        integer, allocatable :: surface_indices(:)
        integer :: n_surfaces = 0
    end type km_d11_result_t

contains

subroutine calc_km_d11_profile(surface_indices, result, &
    per_surface_energy_eV, per_surface_density, per_surface_temperature)

    use constants, only: echarge, ev2erg
    use gorilla_applets_types_mod, only: c, dc, in, s, start
    use km_benchmark_diagnostics_mod, only: collision_trace_t, &
        compute_collision_trace, compute_tracing_time_seconds, &
        estimate_collision_time_seconds
    use km_benchmark_settings_mod, only: boole_run_energy_scan, &
        boole_run_trace_scan, boole_write_surface_trace, diagnostics_prefix, &
        energy_scan_max_factor, energy_scan_min_factor, energy_scan_points, &
        fit_end_fraction, fit_start_fraction, n_trace_scan_multipliers, &
        temperature_eV, trace_scan_multipliers, trace_time_multiplier
    use tetra_grid_mod, only: verts_sthetaphi
    use tetra_grid_settings_mod, only: grid_size
    use transport_benchmark_utils_mod, only: set_active_species_parameters
    use utils_data_pre_and_post_processing_mod, only: &
        calc_collision_coefficients_for_all_tetrahedra, initialize_exit_data
    use utils_self_consistent_ef_mod, only: allocate_start_type, &
        set_particle_type_specifications

    integer, intent(in) :: surface_indices(:)
    type(km_d11_result_t), intent(out) :: result
    real(dp), intent(in), optional :: per_surface_energy_eV(:)
    real(dp), intent(in), optional :: per_surface_density(:, :)
    real(dp), intent(in), optional :: per_surface_temperature(:, :)

    real(dp), allocatable :: background_charge(:)
    real(dp), allocatable :: background_density(:)
    real(dp), allocatable :: background_mass(:)
    real(dp), allocatable :: background_temperature(:)
    real(dp), allocatable :: boundary_s(:)
    real(dp) :: surface_energy
    integer :: i, n_bg, n_particles, ns, idx
    type(collision_trace_t) :: collision_trace

    n_particles = in%num_particles
    if (n_particles < 100) n_particles = 5000

    result%n_surfaces = size(surface_indices)
    allocate(result%surface_indices(result%n_surfaces))
    allocate(result%d11(result%n_surfaces))
    allocate(result%convection(result%n_surfaces))
    allocate(result%energy_eV(result%n_surfaces))
    allocate(result%collision_time_s(result%n_surfaces))
    allocate(result%trace_time_s(result%n_surfaces))
    allocate(result%total_collision_frequency_hz(result%n_surfaces))
    allocate(result%convolved_d11(result%n_surfaces))
    result%surface_indices = surface_indices
    result%d11 = 0.0_dp
    result%convection = 0.0_dp
    result%energy_eV = 0.0_dp
    result%collision_time_s = 0.0_dp
    result%trace_time_s = 0.0_dp
    result%total_collision_frequency_hz = 0.0_dp
    result%convolved_d11 = 0.0_dp

    allocate(boundary_s(grid_size(1) + 1))
    do i = 1, grid_size(1) + 1
        boundary_s(i) = verts_sthetaphi(1, grid_size(3) * (i - 1) + 1)
    end do
    result%boundary_s = boundary_s

    allocate(background_density(size(in%background_density_cm3)))
    allocate(background_temperature(size(in%background_temperature_eV)))
    allocate(background_mass(size(in%background_mass)))
    allocate(background_charge(size(in%background_charge_num)))
    background_density = in%background_density_cm3
    background_temperature = in%background_temperature_eV
    background_mass = in%background_mass
    background_charge = in%background_charge_num
    n_bg = size(background_mass)

    s%temperature = in%energy_eV
    s%n_particles = n_particles
    s%k = 1000
    s%j = 100
    call reset_km_buffers(n_particles)

    if (.not. allocated(dc%s_vertices)) allocate(dc%s_vertices(grid_size(1) + 1))
    dc%s_vertices = boundary_s

    call allocate_start_type(n_particles)
    call set_particle_type_specifications()
    call initialize_exit_data(n_particles)

    do idx = 1, result%n_surfaces
        ns = surface_indices(idx)
        if (ns < 2 .or. ns > grid_size(1)) cycle

        surface_energy = in%energy_eV
        if (present(per_surface_energy_eV)) surface_energy = per_surface_energy_eV(idx)
        if (present(per_surface_density)) then
            background_density = per_surface_density(idx, 1:n_bg)
        end if
        if (present(per_surface_temperature)) then
            background_temperature = per_surface_temperature(idx, 1:n_bg)
        end if

        call set_active_species_parameters(in%tracer_species)
        call apply_surface_collision_state(surface_energy, background_density, &
            background_temperature, background_mass, background_charge)
        call compute_collision_trace(start%particle_mass(in%tracer_species), &
            start%particle_charge(in%tracer_species), background_mass, &
            background_charge, background_density, background_temperature, &
            surface_energy, collision_trace)

        result%energy_eV(idx) = surface_energy
        result%collision_time_s(idx) = estimate_collision_time_seconds(collision_trace)
        result%trace_time_s(idx) = compute_tracing_time_seconds(in%time_step, &
            trace_time_multiplier, result%collision_time_s(idx))
        result%total_collision_frequency_hz(idx) = &
            collision_trace%total_collision_frequency_hz

        print '("  KM D11: surface ",I0,"/",I0," at s=",F8.5," E=",F8.1, &
            &" eV tau_c=",ES12.4," trace=",ES12.4)', idx, result%n_surfaces, &
            boundary_s(ns), surface_energy, result%collision_time_s(idx), &
            result%trace_time_s(idx)

        call run_km_surface(boundary_s(ns), surface_energy, result%trace_time_s(idx), &
            fit_start_fraction, fit_end_fraction, n_particles, result%convection(idx), &
            result%d11(idx))

        if (boole_write_surface_trace) then
            call write_surface_moment_trace(diagnostics_prefix, idx, boundary_s(ns), &
                result%trace_time_s(idx), result%d11(idx), result%convection(idx))
            call write_collision_trace(diagnostics_prefix, idx, boundary_s(ns), &
                background_density, background_temperature, background_charge, &
                collision_trace)
        end if

        result%convolved_d11(idx) = result%d11(idx)
        if (boole_run_trace_scan) then
            call write_trace_scan(diagnostics_prefix, idx, boundary_s(ns), &
                surface_energy, background_density, background_temperature, &
                background_mass, background_charge, n_particles)
        end if
        if (boole_run_energy_scan) then
            call write_energy_scan(diagnostics_prefix, idx, boundary_s(ns), &
                surface_energy, background_density, background_temperature, &
                background_mass, background_charge, n_particles, &
                result%convolved_d11(idx))
        end if

        print '(A,ES12.4,A,ES12.4)', '    D11=', result%d11(idx), '  A=', &
            result%convection(idx)
    end do

contains

    subroutine apply_surface_collision_state(surface_energy, density, temp, mass, charge)

        use utils_data_pre_and_post_processing_mod, only: set_custom_background

        real(dp), intent(in) :: surface_energy
        real(dp), intent(in) :: density(:)
        real(dp), intent(in) :: temp(:)
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: charge(:)

        in%energy_eV = surface_energy
        s%temperature = surface_energy
        start%v0(in%tracer_species) = sqrt(2.0_dp * surface_energy * ev2erg / &
            start%particle_mass(in%tracer_species))
        call set_custom_background(size(mass), density, temp, mass, charge)
        call calc_collision_coefficients_for_all_tetrahedra(in%tracer_species)

    end subroutine apply_surface_collision_state

    subroutine reset_km_buffers(local_n_particles)

        integer, intent(in) :: local_n_particles

        if (allocated(s%delta_s)) deallocate(s%delta_s)
        if (allocated(s%delta_s_squared)) deallocate(s%delta_s_squared)
        if (allocated(s%time)) deallocate(s%time)
        if (allocated(s%f_v)) deallocate(s%f_v)
        if (allocated(s%check)) deallocate(s%check)
        if (allocated(s%boole_large_distance)) deallocate(s%boole_large_distance)
        allocate(s%delta_s(s%k))
        allocate(s%delta_s_squared(s%k))
        allocate(s%time(s%k))
        allocate(s%f_v(s%k, s%j))
        allocate(s%check(s%k))
        allocate(s%boole_large_distance(local_n_particles))

    end subroutine reset_km_buffers

end subroutine calc_km_d11_profile

subroutine run_km_surface(surface_s, surface_energy, trace_time_s, &
    fit_start_fraction, fit_end_fraction, n_particles, convection_a, d11_s)

    use gorilla_applets_types_mod, only: in, s
    use transport_benchmark_utils_mod, only: fit_transport_coefficients
    use utils_data_pre_and_post_processing_mod, only: &
        prepare_next_round_of_parallelised_particle_pushing
    use utils_self_consistent_ef_mod, only: parallelised_particle_pushing, &
        set_rest_of_individual_particle_specifications, set_starting_positions

    real(dp), intent(in) :: surface_s
    real(dp), intent(in) :: surface_energy
    real(dp), intent(in) :: trace_time_s
    real(dp), intent(in) :: fit_start_fraction
    real(dp), intent(in) :: fit_end_fraction
    integer, intent(in) :: n_particles
    real(dp), intent(out) :: convection_a
    real(dp), intent(out) :: d11_s

    real(dp), allocatable :: rand_matrix(:, :, :)
    integer :: i

    allocate(rand_matrix(5, n_particles, 1))

    s%s0 = surface_s
    s%temperature = surface_energy
    call set_trace_time_for_species(trace_time_s, in%tracer_species)
    s%time = [(trace_time_s / real(s%k, dp) * real(i, dp), i = 1, s%k)]
    s%delta_s = 0.0_dp
    s%delta_s_squared = 0.0_dp
    s%f_v = 0.0_dp
    s%check = 0
    s%boole_large_distance = .false.

    call random_number(rand_matrix)
    call set_starting_positions(rand_matrix, (/in%tracer_species/), s%s0)
    call set_rest_of_individual_particle_specifications(rand_matrix, &
        boole_diffusion_coefficient_in=.false., species_in=(/in%tracer_species/), &
        n_particles_in=n_particles)

    call prepare_next_round_of_parallelised_particle_pushing(in%tracer_species)
    call parallelised_particle_pushing(species=in%tracer_species, j=1, &
        boole_diffusion_coefficient=.true., n_particles_in=n_particles)

    call fit_transport_coefficients(s%time, s%delta_s, s%delta_s_squared, &
        convection_a, d11_s, fit_start_fraction=fit_start_fraction, &
        fit_end_fraction=fit_end_fraction)

end subroutine run_km_surface

subroutine set_trace_time_for_species(trace_time_s, species)

    use gorilla_applets_types_mod, only: start

    real(dp), intent(in) :: trace_time_s
    integer, intent(in) :: species

    if (.not. allocated(start%t)) then
        error stop 'set_trace_time_for_species: start%t is not allocated'
    end if
    if (species < 1 .or. species > size(start%t)) then
        error stop 'set_trace_time_for_species: invalid species index'
    end if

    start%t(species) = trace_time_s

end subroutine set_trace_time_for_species

subroutine write_surface_moment_trace(prefix, surface_index, surface_s, &
    trace_time_s, d11_s, convection_a)

    use gorilla_applets_types_mod, only: s

    character(len=*), intent(in) :: prefix
    integer, intent(in) :: surface_index
    real(dp), intent(in) :: surface_s
    real(dp), intent(in) :: trace_time_s
    real(dp), intent(in) :: d11_s
    real(dp), intent(in) :: convection_a

    character(len=512) :: filename
    integer :: file_unit, i

    write(filename, '(A,"_surface_",I3.3,"_moments.csv")') trim(prefix), surface_index
    open(newunit=file_unit, file=trim(filename), status='replace', action='write')
    write(file_unit, '(A)') &
        'time_s,delta_s,delta_s_squared,check,surface_s,trace_time_s,d11_s,convection_a'
    do i = 1, size(s%time)
        write(file_unit, '(8(ES24.16,:,","))') s%time(i), s%delta_s(i), &
            s%delta_s_squared(i), real(s%check(i), dp), surface_s, trace_time_s, d11_s, &
            convection_a
    end do
    close(file_unit)

end subroutine write_surface_moment_trace

subroutine write_collision_trace(prefix, surface_index, surface_s, density, temp, &
    charge, collision_trace)

    use km_benchmark_diagnostics_mod, only: collision_trace_t

    character(len=*), intent(in) :: prefix
    integer, intent(in) :: surface_index
    real(dp), intent(in) :: surface_s
    real(dp), intent(in) :: density(:)
    real(dp), intent(in) :: temp(:)
    real(dp), intent(in) :: charge(:)
    type(collision_trace_t), intent(in) :: collision_trace

    character(len=512) :: filename
    integer :: file_unit, i

    write(filename, '(A,"_surface_",I3.3,"_collision.csv")') trim(prefix), surface_index
    open(newunit=file_unit, file=trim(filename), status='replace', action='write')
    write(file_unit, '(A)') &
        'species_index,surface_s,charge,density_cm3,temperature_eV,lambda_ab,velrat,enrat,efcolf,collision_frequency_hz'
    do i = 1, size(density)
        write(file_unit, '(I0,",",9(ES24.16,:,","))') i, surface_s, charge(i), &
            density(i), temp(i), collision_trace%lambda(i), &
            collision_trace%velrat(i), collision_trace%enrat(i), &
            collision_trace%efcolf(i), collision_trace%collision_frequency_hz(i)
    end do
    close(file_unit)

end subroutine write_collision_trace

subroutine write_trace_scan(prefix, surface_index, surface_s, surface_energy, &
    density, temp, mass, charge, n_particles)

    use gorilla_applets_types_mod, only: in, start
    use km_benchmark_diagnostics_mod, only: collision_trace_t, &
        compute_collision_trace, compute_relative_change, compute_tracing_time_seconds, &
        estimate_collision_time_seconds
    use km_benchmark_settings_mod, only: fit_end_fraction, fit_start_fraction, &
        n_trace_scan_multipliers, trace_scan_multipliers

    character(len=*), intent(in) :: prefix
    integer, intent(in) :: surface_index
    real(dp), intent(in) :: surface_s
    real(dp), intent(in) :: surface_energy
    real(dp), intent(in) :: density(:)
    real(dp), intent(in) :: temp(:)
    real(dp), intent(in) :: mass(:)
    real(dp), intent(in) :: charge(:)
    integer, intent(in) :: n_particles

    type(collision_trace_t) :: collision_trace
    character(len=512) :: filename
    real(dp) :: collision_time_s, d11_s, convection_a, relative_change, trace_time_s
    real(dp) :: previous_d11
    integer :: file_unit, i

    write(filename, '(A,"_surface_",I3.3,"_trace_scan.csv")') trim(prefix), surface_index
    open(newunit=file_unit, file=trim(filename), status='replace', action='write')
    write(file_unit, '(A)') &
        'surface_s,multiplier,collision_time_s,trace_time_s,d11_s,convection_a,relative_change'

    previous_d11 = -1.0_dp
    do i = 1, n_trace_scan_multipliers
        call set_active_surface_state(surface_energy, density, temp, mass, charge)
        call compute_collision_trace(start%particle_mass(in%tracer_species), &
            start%particle_charge(in%tracer_species), mass, charge, density, temp, &
            surface_energy, collision_trace)
        collision_time_s = estimate_collision_time_seconds(collision_trace)
        trace_time_s = compute_tracing_time_seconds(0.0_dp, trace_scan_multipliers(i), &
            collision_time_s)
        call run_km_surface(surface_s, surface_energy, trace_time_s, fit_start_fraction, &
            fit_end_fraction, n_particles, convection_a, d11_s)
        relative_change = 0.0_dp
        if (previous_d11 > 0.0_dp) then
            relative_change = compute_relative_change(previous_d11, d11_s)
        end if
        write(file_unit, '(7(ES24.16,:,","))') surface_s, trace_scan_multipliers(i), &
            collision_time_s, trace_time_s, d11_s, convection_a, relative_change
        previous_d11 = d11_s
    end do

    close(file_unit)

contains

    subroutine set_active_surface_state(active_energy, active_density, active_temp, &
        active_mass, active_charge)

        use constants, only: ev2erg
        use gorilla_applets_types_mod, only: in, s, start
        use utils_data_pre_and_post_processing_mod, only: &
            calc_collision_coefficients_for_all_tetrahedra, set_custom_background

        real(dp), intent(in) :: active_energy
        real(dp), intent(in) :: active_density(:)
        real(dp), intent(in) :: active_temp(:)
        real(dp), intent(in) :: active_mass(:)
        real(dp), intent(in) :: active_charge(:)

        in%energy_eV = active_energy
        s%temperature = active_energy
        start%v0(in%tracer_species) = sqrt(2.0_dp * active_energy * ev2erg / &
            start%particle_mass(in%tracer_species))
        call set_custom_background(size(active_mass), active_density, active_temp, &
            active_mass, active_charge)
        call calc_collision_coefficients_for_all_tetrahedra(in%tracer_species)

    end subroutine set_active_surface_state

end subroutine write_trace_scan

subroutine write_energy_scan(prefix, surface_index, surface_s, surface_energy, &
    density, temp, mass, charge, n_particles, convolved_d11)

    use gorilla_applets_types_mod, only: in, start
    use km_benchmark_diagnostics_mod, only: build_log_energy_grid, &
        collision_trace_t, compute_collision_trace, compute_energy_convolution, &
        compute_tracing_time_seconds, estimate_collision_time_seconds
    use km_benchmark_settings_mod, only: energy_scan_max_factor, &
        energy_scan_min_factor, energy_scan_points, fit_end_fraction, &
        fit_start_fraction, trace_time_multiplier

    character(len=*), intent(in) :: prefix
    integer, intent(in) :: surface_index
    real(dp), intent(in) :: surface_s
    real(dp), intent(in) :: surface_energy
    real(dp), intent(in) :: density(:)
    real(dp), intent(in) :: temp(:)
    real(dp), intent(in) :: mass(:)
    real(dp), intent(in) :: charge(:)
    integer, intent(in) :: n_particles
    real(dp), intent(out) :: convolved_d11

    type(collision_trace_t) :: collision_trace
    character(len=512) :: filename
    real(dp), allocatable :: collision_time_s(:)
    real(dp), allocatable :: convection(:)
    real(dp), allocatable :: d11_s(:)
    real(dp), allocatable :: energy_grid(:)
    real(dp), allocatable :: trace_time_s(:)
    real(dp), allocatable :: weights(:)
    integer :: file_unit, i

    call build_log_energy_grid(surface_energy, energy_scan_min_factor, &
        energy_scan_max_factor, energy_scan_points, energy_grid)
    allocate(collision_time_s(size(energy_grid)))
    allocate(convection(size(energy_grid)))
    allocate(d11_s(size(energy_grid)))
    allocate(trace_time_s(size(energy_grid)))

    write(filename, '(A,"_surface_",I3.3,"_energy_scan.csv")') trim(prefix), surface_index
    open(newunit=file_unit, file=trim(filename), status='replace', action='write')
    write(file_unit, '(A)') &
        'surface_s,energy_eV,collision_time_s,trace_time_s,d11_s,convection_a,weight'

    do i = 1, size(energy_grid)
        call set_active_surface_state(energy_grid(i), density, temp, mass, charge)
        call compute_collision_trace(start%particle_mass(in%tracer_species), &
            start%particle_charge(in%tracer_species), mass, charge, density, temp, &
            energy_grid(i), collision_trace)
        collision_time_s(i) = estimate_collision_time_seconds(collision_trace)
        trace_time_s(i) = compute_tracing_time_seconds(0.0_dp, trace_time_multiplier, &
            collision_time_s(i))
        call run_km_surface(surface_s, energy_grid(i), trace_time_s(i), &
            fit_start_fraction, fit_end_fraction, n_particles, convection(i), d11_s(i))
    end do

    call compute_energy_convolution(energy_grid, d11_s, surface_energy, &
        convolved_d11, weights)
    do i = 1, size(energy_grid)
        write(file_unit, '(7(ES24.16,:,","))') surface_s, energy_grid(i), &
            collision_time_s(i), trace_time_s(i), d11_s(i), convection(i), weights(i)
    end do
    close(file_unit)

contains

    subroutine set_active_surface_state(active_energy, active_density, active_temp, &
        active_mass, active_charge)

        use constants, only: ev2erg
        use gorilla_applets_types_mod, only: in, s, start
        use utils_data_pre_and_post_processing_mod, only: &
            calc_collision_coefficients_for_all_tetrahedra, set_custom_background

        real(dp), intent(in) :: active_energy
        real(dp), intent(in) :: active_density(:)
        real(dp), intent(in) :: active_temp(:)
        real(dp), intent(in) :: active_mass(:)
        real(dp), intent(in) :: active_charge(:)

        in%energy_eV = active_energy
        s%temperature = active_energy
        start%v0(in%tracer_species) = sqrt(2.0_dp * active_energy * ev2erg / &
            start%particle_mass(in%tracer_species))
        call set_custom_background(size(active_mass), active_density, active_temp, &
            active_mass, active_charge)
        call calc_collision_coefficients_for_all_tetrahedra(in%tracer_species)

    end subroutine set_active_surface_state

end subroutine write_energy_scan

subroutine write_km_csv(result, filename)

    type(km_d11_result_t), intent(in) :: result
    character(len=*), intent(in) :: filename

    integer :: csv_unit, i

    open(newunit=csv_unit, file=filename, status='replace', action='write')
    write(csv_unit, '(A)') &
        's,d11_s,convection_A,energy_eV,collision_time_s,trace_time_s,total_collision_frequency_hz,d11_convolved_s'
    do i = 1, result%n_surfaces
        write(csv_unit, '(8(ES24.16,:,","))') &
            result%boundary_s(result%surface_indices(i)), result%d11(i), &
            result%convection(i), result%energy_eV(i), result%collision_time_s(i), &
            result%trace_time_s(i), result%total_collision_frequency_hz(i), &
            result%convolved_d11(i)
    end do
    close(csv_unit)

end subroutine write_km_csv

end module kramers_moyal_transport_mod
