module kramers_moyal_transport_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    public :: calc_km_d11_profile
    public :: km_d11_result_t

    type km_d11_result_t
        real(dp), allocatable :: boundary_s(:)
        real(dp), allocatable :: d11(:)
        real(dp), allocatable :: convection(:)
        integer, allocatable :: surface_indices(:)
        integer :: n_surfaces = 0
    end type km_d11_result_t

contains

subroutine calc_km_d11_profile(surface_indices, result)

    use gorilla_applets_types_mod, only: in, s, dc, start, g, exit_data, ep
    use tetra_grid_settings_mod, only: grid_size, sfc_s_min
    use tetra_grid_mod, only: verts_sthetaphi
    !
    use transport_benchmark_utils_mod, only: fit_transport_coefficients, &
        set_active_species_parameters
    use utils_data_pre_and_post_processing_mod, only: &
        prepare_next_round_of_parallelised_particle_pushing, &
        calc_collision_coefficients_for_all_tetrahedra, initialize_exit_data
    use utils_self_consistent_ef_mod, only: allocate_start_type, &
        parallelised_particle_pushing, set_particle_type_specifications, &
        set_rest_of_individual_particle_specifications, set_starting_positions

    integer, intent(in) :: surface_indices(:)
    type(km_d11_result_t), intent(out) :: result

    real(dp), allocatable :: rand_matrix(:, :, :)
    real(dp), allocatable :: boundary_s(:)
    real(dp) :: A, B, tau_c_ei
    integer :: ns, i, n_particles, idx
    logical :: first_surface

    n_particles = in%num_particles
    if (n_particles < 100) n_particles = 5000

    result%n_surfaces = size(surface_indices)
    allocate(result%surface_indices(result%n_surfaces))
    allocate(result%d11(result%n_surfaces))
    allocate(result%convection(result%n_surfaces))
    result%surface_indices = surface_indices

    ! Get boundary s-coordinates
    allocate(boundary_s(grid_size(1) + 1))
    do i = 1, grid_size(1) + 1
        boundary_s(i) = verts_sthetaphi(1, grid_size(3) * (i - 1) + 1)
    end do
    result%boundary_s = boundary_s

    ! Use tracing time from input (time_step), or estimate from collision time
    if (in%time_step > 0.0_dp) then
        tau_c_ei = in%time_step
    else
        tau_c_ei = 1.7d-4
    end if

    s%temperature = in%energy_eV

    ! Set up displacement accumulation
    s%n_particles = n_particles
    s%k = 1000
    s%j = 100
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
    allocate(s%boole_large_distance(n_particles))

    ! Set up dc structure for vertex coordinates
    if (.not. allocated(dc%s_vertices)) allocate(dc%s_vertices(grid_size(1) + 1))
    dc%s_vertices = boundary_s

    allocate(rand_matrix(5, n_particles, 1))

    call allocate_start_type(n_particles)
    call set_particle_type_specifications()
    call initialize_exit_data(n_particles)

    ! Tracing time: 2x collision time (must be after allocate_start_type)
    start%t(in%tracer_species) = 2.0_dp * tau_c_ei
    s%time = [(start%t(in%tracer_species) / s%k * i, i = 1, s%k)]

    first_surface = .true.

    do idx = 1, result%n_surfaces
        ns = surface_indices(idx)
        if (ns < 2 .or. ns > grid_size(1)) cycle

        print '(A,I0,A,I0,A,F8.5)', '  KM D11: surface ', idx, '/', &
            result%n_surfaces, ' at s=', boundary_s(ns)

        s%s0 = boundary_s(ns)
        s%delta_s = 0.0_dp
        s%delta_s_squared = 0.0_dp
        s%f_v = 0.0_dp
        s%check = 0
        s%boole_large_distance = .false.

        call random_number(rand_matrix)
        call set_starting_positions(rand_matrix, (/in%tracer_species/), s%s0)
        call set_rest_of_individual_particle_specifications(rand_matrix, &
            boole_diffusion_coefficient_in=.false., &
            species_in=(/in%tracer_species/), n_particles_in=n_particles)

        call prepare_next_round_of_parallelised_particle_pushing(in%tracer_species)

        if (first_surface) then
            call calc_collision_coefficients_for_all_tetrahedra(in%tracer_species)
            first_surface = .false.
        end if

        call parallelised_particle_pushing(species=in%tracer_species, j=1, &
            boole_diffusion_coefficient=.true., n_particles_in=n_particles)

        call fit_transport_coefficients(s%time, s%delta_s, s%delta_s_squared, A, B)

        result%d11(idx) = B
        result%convection(idx) = A

        print '(A,ES12.4,A,ES12.4)', '    D11=', B, '  A=', A
    end do

end subroutine calc_km_d11_profile

end module kramers_moyal_transport_mod
