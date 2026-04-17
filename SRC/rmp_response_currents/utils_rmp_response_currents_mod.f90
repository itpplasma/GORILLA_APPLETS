module utils_rmp_response_currents_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine read_rmp_response_currents_inp_into_type

    use gorilla_applets_types_mod, only: in

    real(dp) :: time_step, energy_eV, n_particles, density
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_monoenergetic, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, &
               boole_write_exit_data, boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, &
               boole_eliminate_particles_outside_flux, boole_delta_f
    integer :: i_integrator_type, seed_option, n_species
    real(dp) :: flux_threshold_for_elimination

    integer :: s_inp_unit

    NAMELIST /rmp_response_currents_nml/ time_step, energy_eV, n_particles, boole_squared_moments, boole_point_source, &
    & boole_collisions, boole_precalc_collisions, density, boole_refined_sqrt_g, boole_monoenergetic, &
    & boole_linear_density_simulation, boole_antithetic_variate, boole_linear_temperature_simulation, i_integrator_type, &
    & seed_option, boole_write_vertex_indices, boole_write_vertex_coordinates, boole_write_prism_volumes, &
    & boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_species, &
    & boole_eliminate_particles_outside_flux, flux_threshold_for_elimination, boole_delta_f

    open(newunit = s_inp_unit, file='rmp_response_currents.inp', status='unknown')
    read(s_inp_unit,nml=rmp_response_currents_nml)
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
    in%boole_monoenergetic = boole_monoenergetic
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
    in%boole_eliminate_particles_outside_flux = boole_eliminate_particles_outside_flux
    in%flux_threshold_for_elimination = flux_threshold_for_elimination
    in%boole_delta_f = boole_delta_f

end subroutine read_rmp_response_currents_inp_into_type

! ====================================================================
subroutine parallelised_particle_pushing_rmp_response_currents(species, n_particles_in)

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

            if (start%lost(n, species)) then
                call update_exit_data(.true., 0.0_dp, start%x(:,n,species), 0.0_dp, 0.0_dp, 0, n, species_in=species, ind_tetr=-1)
                cycle
            endif

            call initialise_loop_variables(l, n, local_counter, particle_status, t, local_tetr_moments, x, vpar, vperp, species)

            i = 0

            do while (t%confined.lt.start%t(species))
                i = i + 1

                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, t, x, vpar, vperp, ind_tetr, iface, species, iswmode_in=1)
                    t%step = t%step / start%v0(species)
                else
                    t%step = start%t(species) - t%confined
                endif

                call orbit_timestep_rmp_response_currents(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                                          local_tetr_moments, local_counter, species, start%t(species))

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

end subroutine parallelised_particle_pushing_rmp_response_currents

! ====================================================================
subroutine orbit_timestep_rmp_response_currents(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                                local_tetr_moments, local_counter, species, t_tot)

    use pusher_tetra_rk_mod, only: pusher_tetra_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: vperp_func
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_types_mod, only: counter_t, particle_status_t, start, in, time_t, g, weights
    use tetra_grid_settings_mod, only: grid_kind, sfc_s_min
    use utils_orbit_timestep_mod, only: initialize_constants_of_motion, compute_radial_fluxes, &
        identify_particles_entering_annulus, update_local_tetr_moments
    use utils_helical_core_mod, only: apply_weight_fading

    real(dp), dimension(3), intent(inout)        :: x
    real(dp), intent(inout)                      :: vpar, vperp
    type(time_t), intent(inout)                  :: t
    type(particle_status_t), intent(inout)       :: particle_status
    integer, intent(inout)                       :: ind_tetr, iface
    integer, intent(in)                          :: n, species
    complex(dp), dimension(:,:), intent(inout)   :: local_tetr_moments
    type(counter_t), intent(inout)               :: local_counter
    real(dp), intent(in)                         :: t_tot

    real(dp), dimension(3)                       :: z_save, x_new
    real(dp)                                     :: t_pass, perpinv, rand_frac
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
        call calc_particle_weights_and_jperp_rmp_response_currents(n, z_save, vpar, vperp, ind_tetr, species)
        if (in%boole_delta_f) weights%original(n, species) = weights%w(n, species)
        particle_status%initialized = .true.
    endif

    if (t%step.eq.0.0_dp) return
    if (particle_status%initialized) z_save = x - tetra_physics(ind_tetr)%x1
    call initialize_constants_of_motion(vperp, z_save, ind_tetr, perpinv)

    t%remain = t%step
    boole_t_finished = .false.
    local_counter%tetr_pushings = local_counter%tetr_pushings - 1

    do
        local_counter%tetr_pushings = local_counter%tetr_pushings + 1

        if (ind_tetr.eq.-1) then
            if ((grid_kind.eq.2).or.(grid_kind.eq.3)) then
                call identify_particles_entering_annulus(x, local_counter, boole_lost_inside)
                if (boole_lost_inside) then
                    x_new = 3*(/g%raxis, x(2), g%zaxis/) - 2*x
                    vperp = vperp_func(z_save, perpinv, ind_tetr_save)
                    call find_tetra(x_new, vpar, vperp, ind_tetr, iface)
                    x = x_new
                    if (ind_tetr.eq.-1) then
                        print*, "ATTENTION: particle pushing across the hole surrounding the magnetic axis was unsuccessful"
                        exit
                    endif
                else
                    call random_number(rand_frac)
                    rand_frac = 0.01_dp * rand_frac
                    x_new(1) = x(1) + rand_frac * (g%raxis - x(1))
                    x_new(2) = x(2)
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

        if (in%boole_delta_f) call apply_weight_fading(n, species, t, t_tot)

        call update_local_tetr_moments(local_tetr_moments, ind_tetr_save, n, optional_quantities, species)
        if ((grid_kind.eq.2).or.(grid_kind.eq.3)) call compute_radial_fluxes(ind_tetr_save, ind_tetr, x)

        if (boole_t_finished) exit
    enddo

end subroutine orbit_timestep_rmp_response_currents

! ====================================================================
subroutine calc_particle_weights_and_jperp_rmp_response_currents(n, z_save, vpar, vperp, ind_tetr, species_in)

    use gorilla_applets_types_mod, only: in, flux, start, weights, g
    use tetra_physics_mod, only: tetra_physics
    use constants, only: ev2erg, pi
    use volume_integrals_and_sqrt_g_mod, only: sqrt_g
    use supporting_functions_mod, only: bmod_func
    use gorilla_settings_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind
    use marker_distribution_mod, only: evaluate_distribution_3d, evaluate_distribution_1d, pdf_boltzmann
    use utils_helical_core_mod, only: adapt_weights_delta_f

    real(dp), intent(in) :: vpar, vperp
    real(dp), dimension(3), intent(in) :: z_save
    integer, intent(in) :: n, ind_tetr
    integer, intent(in), optional :: species_in
    integer :: species = 1
    real(dp) :: local_poloidal_flux
    real(dp) :: r, phi, z
    real(dp) :: s_value
    real(dp) :: f, J_x, J_y
    real(dp) :: x_global(3), pdf_position, pdf_energy_ev, pdf_energy_erg, pdf_lambda, d_lambda_epsilon_d_jperp_vpar, omega_c, &
    pdf_velocity_space, pdf_gyroangle
    real(dp) :: density, T, m, energy_ev, energy_erg, q, Phi_elec

    if (present(species_in)) species = species_in

    r = z_save(1)
    phi = z_save(2)
    z = z_save(3)

    if (in%boole_linear_density_simulation.or.in%boole_linear_temperature_simulation) then
        if (coord_system == 2) then
            s_value = tetra_physics(ind_tetr)%x1(1) + z_save(1)
        else if (grid_kind /= 3) then
            local_poloidal_flux = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*z_save)
            s_value = local_poloidal_flux / flux%poloidal_max
        else
            print*, 'Error in calc_particle_weights_and_jperp_rmp_response_currents: Computing radial coordinate from A_phi is &
                    &only valid for axisymmetric devices. For stellarators (grid_kind=3), use flux coordinates (coord_system=2).'
            stop
        endif
    endif

    m = start%particle_mass(species)
    T = in%energy_eV*ev2erg
    energy_ev = start%energy(n,species)
    energy_erg = start%energy(n,species)*ev2erg
    density = in%density

    x_global = z_save + tetra_physics(ind_tetr)%x1
    pdf_position = evaluate_distribution_3d(start%dist_position, x_global)
    pdf_energy_ev = evaluate_distribution_1d(start%dist_energy, energy_ev)
    if (pdf_energy_ev.eq.1.0_dp) then
        pdf_energy_erg = 1.0_dp
    else
        pdf_energy_erg = pdf_energy_ev/ev2erg
    endif
    pdf_lambda = evaluate_distribution_1d(start%dist_lambda, start%pitch(n,species))

    if (in%boole_refined_sqrt_g) then
        J_x = (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
            & (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))
    else
        J_x = r + tetra_physics(ind_tetr)%x1(1)
    endif

    if ((.not.in%boole_delta_f).and.(in%boole_linear_density_simulation)) then
        density = density*(1.1_dp - s_value)/1.1_dp
    endif

    pdf_gyroangle = 1/(2*pi)

    if (in%boole_monoenergetic) then
        f = density*pdf_gyroangle/(2.0_dp*start%v0(species))
        J_y = start%v0(species)
        pdf_velocity_space = pdf_lambda*pdf_gyroangle
    else
        q = start%particle_charge(species)
        Phi_elec = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z_save)
        if (in%boole_linear_temperature_simulation) then
            T = T*(1.1_dp - s_value)/1.1_dp
        endif
        f = pdf_boltzmann(density, T, m, energy_erg, q, Phi_elec)
        omega_c = bmod_func(z_save, ind_tetr)/start%cm_over_e(species)
        d_lambda_epsilon_d_jperp_vpar = omega_c*sqrt(m/(2.0_dp*energy_erg))
        J_y = m**2*omega_c
        pdf_velocity_space = pdf_lambda*pdf_energy_erg*pdf_gyroangle*d_lambda_epsilon_d_jperp_vpar
    endif
    if (in%boole_point_source) f = f*((g%amax-g%amin)*(g%cmax-g%cmin)*2*pi)

    weights%w(n,species) = f*J_x*J_y/(pdf_position*pdf_velocity_space)
    start%jperp(n,species) = m*vperp**2*start%cm_over_e(species)/(2*bmod_func(z_save,ind_tetr))*(-1)

    if (in%boole_delta_f) call adapt_weights_delta_f(n, z_save, vpar, vperp, ind_tetr, species)

end subroutine calc_particle_weights_and_jperp_rmp_response_currents

end module utils_rmp_response_currents_mod
