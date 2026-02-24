module utils_helical_core_mod
!
! Module for helical core particle tracing utilities.
!
! This is a simplified version of utils_anomalous_transport_mod that:
!   - Provides particle pushing without anomalous displacement
!   - Does NOT include straight-line displacement routines
!   - Does NOT include electric potential perturbation scan routines
!   - Does NOT include diffusion coefficient calculation routines
!
! Contains routines for:
!   - Reading helical core input parameters
!   - Parallelised particle pushing (without anomalous displacement)
!   - Orbit timestep integration (standard guiding-center)
!

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine read_helical_core_inp_into_type

    use gorilla_applets_types_mod, only: in

    real(dp) :: time_step, energy_eV, n_particles, density
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_boltzmann_energies, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, &
               boole_write_exit_data, boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions
    integer :: i_integrator_type, seed_option, n_species

    integer :: s_inp_unit

    NAMELIST /helical_core_nml/ time_step, energy_eV, n_particles, boole_squared_moments, boole_point_source, &
    & boole_collisions, boole_precalc_collisions, density, boole_refined_sqrt_g, boole_boltzmann_energies, &
    & boole_linear_density_simulation, boole_antithetic_variate, boole_linear_temperature_simulation, i_integrator_type, &
    & seed_option, boole_write_vertex_indices, boole_write_vertex_coordinates, boole_write_prism_volumes, &
    & boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_species

    open(newunit = s_inp_unit, file='helical_core.inp', status='unknown')
    read(s_inp_unit,nml=helical_core_nml)
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
    in%boole_write_moments = boole_write_moments
    in%boole_write_fourier_moments = boole_write_fourier_moments
    in%boole_write_exit_data = boole_write_exit_data
    in%boole_write_grid_data = boole_write_grid_data
    in%boole_preserve_energy_and_momentum_during_collisions = boole_preserve_energy_and_momentum_during_collisions
    in%n_species = n_species

    ! Set defaults for anomalous transport-specific parameters (disabled)
    in%boole_calc_diffusion_coefficient = .false.
    in%i_scan_option = 0
    in%anomalous_diffusion_coefficient = 0.0_dp

    print *,'GORILLA_APPLETS: Loaded input data from helical_core.inp'

end subroutine read_helical_core_inp_into_type

! ====================================================================
subroutine parallelised_particle_pushing_helical_core(species, n_particles_in)
!
! Parallelised particle pushing for helical core.
! Standard guiding-center orbit integration without anomalous displacement.
!
    use gorilla_applets_types_mod, only: counter, c, in, time_t, moment_specs, counter_t, particle_status_t, start, s
    use tetra_grid_mod, only: ntetr
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use utils_parallelised_particle_pushing_mod, only: print_progress, handle_lost_particles, add_local_tetr_moments_to_output, &
        add_local_counter_to_counter, initialise_loop_variables, carry_out_collisions, update_exit_data, update_start_type, &
        initialise_seed_for_random_numbers_for_each_thread

    integer, intent(in)                               :: species
    integer, intent(in), optional                     :: n_particles_in
    integer                                           :: kpart, iantithetic, ind_tetr, iface, n_particles
    integer                                           :: p, l, n, i
    real(dp), dimension(3)                            :: x
    real(dp)                                          :: vpar, vperp, t_tot
    type(time_t)                                      :: t
    type(counter_t)                                   :: local_counter
    type(particle_status_t)                           :: particle_status
    complex(dp), dimension(:,:), allocatable          :: local_tetr_moments
    logical                                           :: thread_flag = .true.

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

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
                else
                    t%step = start%t(species) - t%confined
                endif

                ! Perform guiding-center orbit integration (NO anomalous displacement)
                call orbit_timestep_helical_core(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
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

            call update_exit_data(particle_status%lost, t%confined, x, vpar, vperp, i, n, species_in=species, ind_tetr=ind_tetr)
            call update_start_type(x, vpar, vperp, n, species, ind_tetr)
        enddo

        !$omp critical
        call add_local_tetr_moments_to_output(local_tetr_moments, species)
        !$omp end critical
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    print*, 'Total tracing time / number of particles: ', t_tot/n_particles, 's'

end subroutine parallelised_particle_pushing_helical_core

! ====================================================================
subroutine orbit_timestep_helical_core(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                    local_tetr_moments, local_counter, species)
!
! Standard orbit timestep for helical core.
! Performs guiding-center orbit integration without anomalous displacement.
!
    use pusher_tetra_rk_mod, only: pusher_tetra_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: vperp_func
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_types_mod, only: counter_t, particle_status_t, start, in, time_t, g
    use tetra_grid_settings_mod, only: grid_kind, sfc_s_min
    use utils_orbit_timestep_mod, only: update_local_tetr_moments, initialize_constants_of_motion, compute_radial_fluxes, &
        calc_particle_weights_and_jperp, identify_particles_entering_annulus

    real(dp), dimension(3), intent(inout)        :: x
    real(dp), intent(inout)                      :: vpar, vperp
    type(time_t), intent(inout)                  :: t
    type(particle_status_t), intent(inout)       :: particle_status
    integer, intent(inout)                       :: ind_tetr, iface
    integer, intent(in)                          :: n, species
    complex(dp), dimension(:,:), intent(inout)   :: local_tetr_moments
    type(counter_t), intent(inout)               :: local_counter

    real(dp), dimension(3)                       :: z_save, x_new
    real(dp)                                     :: t_pass, perpinv
    real(dp)                                     :: rand_frac
    logical                                      :: boole_t_finished, boole_lost_inside
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
        call calc_particle_weights_and_jperp(n, z_save, vpar, vperp, ind_tetr, species)
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
            if ((grid_kind.eq.2).or.(grid_kind.eq.3)) then
                call identify_particles_entering_annulus(x, local_counter, boole_lost_inside)
                if (boole_lost_inside) then
                    ! Particle is near the magnetic axis - reflect across it
                    x_new = 3*(/g%raxis, x(2), g%zaxis/) - 2*x
                    vperp = vperp_func(z_save, perpinv, ind_tetr_save)
                    call find_tetra(x_new, vpar, vperp, ind_tetr, iface)
                    x = x_new
                    if (ind_tetr.eq.-1) then
                        print*, "ATTENTION: particle pushing across the hole surrounding the magnetic axis was unsuccessful"
                        exit
                    else
                        print*, "particle pushing across the hole surrounding the magnetic axis was successful"
                    endif
                else
                    ! Particle left at the outer boundary - displace toward the magnetic axis
                    ! Keep phi unchanged, move R and Z toward axis by random fraction (0 to 1%) of distance
                    call random_number(rand_frac)
                    rand_frac = 0.01_dp * rand_frac  ! 0 to 1% of distance to axis
                    x_new(1) = x(1) + rand_frac * (g%raxis - x(1))
                    x_new(2) = x(2)  ! phi unchanged
                    x_new(3) = x(3) + rand_frac * (g%zaxis - x(3))
                    vperp = vperp_func(z_save, perpinv, ind_tetr_save)
                    call find_tetra(x_new, vpar, vperp, ind_tetr, iface)
                    if (ind_tetr.ne.-1) then
                        x = x_new
                    else
                        exit
                    endif
                endif
            else
                exit
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

end subroutine orbit_timestep_helical_core

end module utils_helical_core_mod
