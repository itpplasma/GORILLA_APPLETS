module utils_anomalous_transport_mod
!
! Module for anomalous transport calculations.
!
! Contains routines for:
!   - Reading anomalous transport input parameters
!   - Parallelised particle pushing (simplified, without self-consistent EF)
!   - Orbit timestep integration (simplified)
!   - Straight-line radial displacement for anomalous diffusion
!

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine read_anomalous_transport_inp_into_type

    use gorilla_applets_types_mod, only: in

    real(dp) :: time_step,energy_eV,n_particles, density
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_boltzmann_energies, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
               boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
               boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, boole_static_ne
    integer :: i_integrator_type, seed_option, n_electric_potential_updates, update_dimension, n_species

    integer :: s_inp_unit

    !Namelist for anomalous transport input (uses same namelist name for compatibility)
    NAMELIST /self_consistent_ef_nml/ time_step,energy_eV,n_particles,boole_squared_moments,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, &
    & boole_antithetic_variate,boole_linear_temperature_simulation,i_integrator_type,seed_option, boole_write_vertex_indices, &
    & boole_write_vertex_coordinates, boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
    & boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_electric_potential_updates, update_dimension, &
    & n_species, boole_static_ne

    open(newunit = s_inp_unit, file='anomalous_transport.inp', status='unknown')
    read(s_inp_unit,nml=self_consistent_ef_nml)
    close(s_inp_unit)

    in%time_step = time_step
    in%energy_eV = energy_eV
    in%n_particles = n_particles
    in%density = density
    in%boole_squared_moments = boole_squared_moments
    in%boole_point_source = boole_point_source
    in%boole_collisions = boole_collisions
    in%boole_precalc_collisions = boole_precalc_collisions
    in%boole_refined_sqrt_g = boole_refined_sqrt_g
    in%boole_boltzmann_energies = boole_boltzmann_energies
    in%boole_linear_density_simulation = boole_linear_density_simulation
    in%boole_antithetic_variate = boole_antithetic_variate
    in%boole_linear_temperature_simulation = boole_linear_temperature_simulation
    in%i_integrator_type = i_integrator_type
    in%seed_option = seed_option
    in%num_particles = int(n_particles)
    in%boole_write_vertex_indices = boole_write_vertex_indices
    in%boole_write_vertex_coordinates = boole_write_vertex_coordinates
    in%boole_write_prism_volumes = boole_write_prism_volumes
    in%boole_write_refined_prism_volumes = boole_write_refined_prism_volumes
    in%boole_write_boltzmann_density = boole_write_boltzmann_density
    in%boole_write_electric_potential = boole_write_electric_potential
    in%boole_write_moments = boole_write_moments
    in%boole_write_fourier_moments = boole_write_fourier_moments
    in%boole_write_exit_data = boole_write_exit_data
    in%boole_write_grid_data = boole_write_grid_data
    in%boole_preserve_energy_and_momentum_during_collisions = boole_preserve_energy_and_momentum_during_collisions
    in%n_electric_potential_updates = n_electric_potential_updates
    in%update_dimension = update_dimension
    in%n_species = n_species
    in%boole_static_ne = boole_static_ne

    print *,'GORILLA_APPLETS: Loaded input data from anomalous_transport.inp'

end subroutine read_anomalous_transport_inp_into_type

! ====================================================================
subroutine parallelised_particle_pushing_anomalous(species)
!
! Simplified parallelised particle pushing for anomalous transport.
! Removes self-consistent electric field specific parts.
! Adds straight-line radial displacement after collisions.
!
    use gorilla_applets_types_mod, only: counter, c, in, time_t, moment_specs, counter_t, particle_status_t, start, s
    use tetra_grid_mod, only: ntetr
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use utils_parallelised_particle_pushing_mod, only: print_progress, handle_lost_particles, add_local_tetr_moments_to_output, &
        add_local_counter_to_counter, initialise_loop_variables, carry_out_collisions, update_exit_data, update_start_type, &
        initialise_seed_for_random_numbers_for_each_thread

    integer, intent(in)                               :: species
    integer                                           :: kpart, iantithetic, ind_tetr, iface, n_particles
    integer                                           :: p, l, n, i
    real(dp), dimension(3)                            :: x
    real(dp)                                          :: vpar, vperp, t_tot
    type(time_t)                                      :: t
    type(counter_t)                                   :: local_counter
    type(particle_status_t)                           :: particle_status
    complex(dp), dimension(:,:), allocatable          :: local_tetr_moments
    logical                                           :: thread_flag = .true.

    n_particles = in%num_particles

    allocate(local_tetr_moments(moment_specs%n_moments,ntetr))
    kpart = 0
    iantithetic = 1
    if (in%boole_antithetic_variate) iantithetic = 2

    t_tot = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(counter, kpart, species, in, c, iantithetic, start, s, n_particles) &
    !$OMP& REDUCTION(+:t_tot) &
    !$OMP& PRIVATE(p, l, n, i, x, vpar, vperp, t, ind_tetr, iface, local_tetr_moments, local_counter, particle_status) &
    !$OMP& FIRSTPRIVATE(thread_flag)

    if (omp_get_thread_num().eq.0) print*, 'Number of threads: ', omp_get_num_threads()

    !$OMP DO SCHEDULE(static)
    do p = 1, n_particles/iantithetic

        if ((.not.in%boole_precalc_collisions).and.thread_flag) then
            call initialise_seed_for_random_numbers_for_each_thread(omp_get_thread_num(), 1)
            thread_flag = .false.
        endif

        do l = 1, iantithetic
            n = (p-1)*iantithetic + l

            !$omp atomic update
            kpart = kpart + 1
            call print_progress(n_particles, kpart, n)

            call initialise_loop_variables(l, n, local_counter, particle_status, t, local_tetr_moments, x, vpar, vperp, species)

            i = 0

            do while (t%confined.lt.start%t(species))
                i = i + 1

                ! Apply collisions if enabled
                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, t, x, vpar, vperp, ind_tetr, iface, species)
                    t%step = t%step / start%v0(species)

                    ! Apply straight-line radial displacement after collision
                    !call displace_by_straight_line(x, ind_tetr, iface, n, species)
                endif

                ! Perform guiding-center orbit integration
                call orbit_timestep_anomalous(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                              local_tetr_moments, local_counter, species)

                t%confined = t%confined + t%step - t%remain
                t_tot = t_tot + t%step - t%remain

                if (ind_tetr.eq.-1) then
                    call handle_lost_particles(local_counter, particle_status%lost)
                    exit
                endif
            enddo

            !$omp critical
            counter%integration_steps = counter%integration_steps + i
            c%maxcol = max(dble(i)/dble(c%randcoli), c%maxcol)
            call add_local_counter_to_counter(local_counter)
            !$omp end critical

            call update_exit_data(particle_status%lost, t%confined, x, vpar, vperp, i, n, species_in=species)
            call update_start_type(x, vpar, vperp, n, species, ind_tetr)
        enddo

        !$omp critical
        call add_local_tetr_moments_to_output(local_tetr_moments, species)
        !$omp end critical
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    print*, 'Total tracing time / number of particles: ', t_tot/n_particles, 's'

end subroutine parallelised_particle_pushing_anomalous

! ====================================================================
subroutine orbit_timestep_anomalous(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                    local_tetr_moments, local_counter, species)
!
! Simplified orbit timestep for anomalous transport.
! Performs guiding-center orbit integration without self-consistent EF specific parts.
!
    use pusher_tetra_rk_mod, only: pusher_tetra_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: vperp_func
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_types_mod, only: counter_t, particle_status_t, start, in, time_t
    use tetra_grid_settings_mod, only: grid_kind, sfc_s_min
    use utils_orbit_timestep_mod, only: update_local_tetr_moments, initialize_constants_of_motion, compute_radial_fluxes
    use utils_self_consistent_ef_mod, only: calc_particle_weights_and_jperp

    real(dp), dimension(3), intent(inout)        :: x
    real(dp), intent(inout)                      :: vpar, vperp
    type(time_t), intent(inout)                  :: t
    type(particle_status_t), intent(inout)       :: particle_status
    integer, intent(inout)                       :: ind_tetr, iface
    integer, intent(in)                          :: n, species
    complex(dp), dimension(:,:), intent(inout)   :: local_tetr_moments
    type(counter_t), intent(inout)               :: local_counter

    real(dp), dimension(3)                       :: z_save
    real(dp)                                     :: t_pass, perpinv
    logical                                      :: boole_t_finished
    integer                                      :: ind_tetr_save, iper_phi
    type(optional_quantities_type)               :: optional_quantities

    if (.not.particle_status%initialized) then
        call check_coordinate_domain(x)
        call find_tetra(x, vpar, vperp, ind_tetr, iface)
        if (ind_tetr.eq.-1) then
            t%remain = t%step
            return
        endif
        z_save = x - tetra_physics(ind_tetr)%x1
        call calc_particle_weights_and_jperp(n, z_save, vpar, vperp, ind_tetr, species, .false.)
        particle_status%initialized = .true.
    endif

    if (t%step.eq.0.0_dp) return
    if (particle_status%initialized) z_save = x - tetra_physics(ind_tetr)%x1
    call initialize_constants_of_motion(vperp, z_save, ind_tetr, perpinv)

    t%remain = t%step
    boole_t_finished = .false.
    local_counter%tetr_pushings = local_counter%tetr_pushings - 1

    do ! Loop for tetrahedron pushings until t%step is reached
        local_counter%tetr_pushings = local_counter%tetr_pushings + 1

        if (ind_tetr.eq.-1) then
            ! Simple boundary handling: reflect at inner boundary, exit at outer
            if (x(1).lt.1.01_dp*sfc_s_min) then
                ! Reflect at inner boundary
                x(1) = 2.0_dp*sfc_s_min - x(1)
                vpar = -vpar
                call find_tetra(x, vpar, vperp, ind_tetr, iface)
                if (ind_tetr.eq.-1) exit
            else
                exit ! Lost at outer boundary
            endif
        endif

        ind_tetr_save = ind_tetr

        select case(ipusher)
            case(1)
                call pusher_tetra_rk(ind_tetr, iface, x, vpar, z_save, t%remain, t_pass, boole_t_finished, iper_phi)
            case(2)
                call pusher_tetra_poly(poly_order, ind_tetr, iface, x, vpar, z_save, t%remain, &
                                       t_pass, boole_t_finished, iper_phi, optional_quantities)
        end select

        vperp = vperp_func(z_save, perpinv, ind_tetr_save)
        t%remain = t%remain - t_pass

        call update_local_tetr_moments(local_tetr_moments, ind_tetr_save, n, optional_quantities, species)
        if ((grid_kind.eq.2).or.(grid_kind.eq.3)) call compute_radial_fluxes(ind_tetr_save, ind_tetr, x)

        if (boole_t_finished) exit
    enddo

end subroutine orbit_timestep_anomalous

! ====================================================================
subroutine displace_by_straight_line(x, ind_tetr, iface, n, species)
!
! Applies a random radial displacement to simulate anomalous transport.
!
! The displacement is purely radial (in the s-coordinate direction),
! with random sign (inward or outward) and magnitude based on
! a diffusion coefficient.
!
! Input/Output:
!   x(3)      - Particle position in (s, theta, phi) coordinates
!   ind_tetr  - Current tetrahedron index (updated after displacement)
!   iface     - Current face index (updated after displacement)
!   n         - Particle index
!   species   - Species index
!
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_types_mod, only: start, in
    use tetra_grid_settings_mod, only: sfc_s_min

    real(dp), dimension(3), intent(inout) :: x
    integer, intent(inout)                :: ind_tetr, iface
    integer, intent(in)                   :: n, species

    real(dp) :: delta_s, rand_sign, rand_uniform
    real(dp) :: vpar_dummy, vperp_dummy
    real(dp), parameter :: diffusion_coefficient = 1.0d-4  ! TODO: make this an input parameter

    ! Generate random radial displacement
    ! delta_s ~ sqrt(2 * D * dt) for diffusive process
    call random_number(rand_uniform)
    rand_sign = sign(1.0_dp, rand_uniform - 0.5_dp)

    call random_number(rand_uniform)
    delta_s = rand_sign * sqrt(2.0_dp * diffusion_coefficient * in%time_step) * sqrt(-2.0_dp * log(rand_uniform))

    ! Apply radial displacement
    x(1) = x(1) + delta_s

    ! Reflect at inner boundary if necessary
    if (x(1).lt.sfc_s_min) then
        x(1) = 2.0_dp*sfc_s_min - x(1)
    endif

    ! Clamp at outer boundary
    if (x(1).gt.1.0_dp) then
        x(1) = 1.0_dp
    endif

    ! Find the new tetrahedron after displacement
    vpar_dummy = 0.0_dp
    vperp_dummy = 0.0_dp
    call find_tetra(x, vpar_dummy, vperp_dummy, ind_tetr, iface)

end subroutine displace_by_straight_line

end module utils_anomalous_transport_mod
