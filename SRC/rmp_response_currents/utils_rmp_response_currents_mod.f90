module utils_rmp_response_currents_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    ! Module-level storage for delta-f weight-evolution inputs.
    ! These are read by read_rmp_response_currents_inp_into_type and
    ! consumed by calc_rmp_response_currents after grid build.
    character(len=512), public :: profile_dir          = './profiles'
    character(len=512), public :: equil_mapping_file   = './flux_functions.dat'
    logical,            public :: boole_constant_delta_B_psi = .true.
    real(dp),           public :: delta_B_psi_const    = 0.0_dp
    integer,            public :: pert_m_fourier       = 0
    integer,            public :: pert_n_fourier       = 0
    character(len=512), public :: delta_B_psi_file     = ''
    integer,            public :: species_for_delta_f  = 1
    ! Diagnostic: if .true., skip the exp(i*(m theta + n phi)) factor in
    ! the constant-amplitude perturbation and return delta_B_psi_const as
    ! a truly uniform (axisymmetric) real-valued perturbation.
    logical,            public :: boole_skip_phase_for_test = .false.

    ! Regularisation parameters (Albert 2016 Eq. 4 with linear damping).
    ! The integration window per particle is N * tau_c long, the damping
    ! switches on at M * tau_c, and nu_r = nu_r_frac * nu_c.
    real(dp),           public :: nu_r_frac                = 0.5_dp
    integer,            public :: n_collision_times_trace  = 10
    integer,            public :: m_collision_times_reg_on = 3
    ! Coulomb logarithm for the local collision frequency.
    real(dp),           public :: coulomb_log              = 17.0_dp

    ! Optional sanity DFT over toroidal mode numbers. With the complex
    ! drive the prism moment is already the (m,n) Fourier amplitude, so
    ! the DFT is only useful as a leakage check.
    logical,            public :: boole_compute_n_modes_dft = .false.

    ! Outer computational-domain cut in flux label s (= psi_tor/psi_tor_edge).
    ! Particles are killed on first crossing of this surface during the trace.
    ! Default 1.0 => no cut, full domain.
    real(dp),           public :: s_outer_cut = 1.0_dp

    ! Bias the starting-position sampler to draw only from the s-window
    ! [s_inner_sample, s_outer_sample]. Default = full plasma.
    real(dp),           public :: s_inner_sample = 0.0_dp
    real(dp),           public :: s_outer_sample = 1.0_dp

    ! Per-particle regularisation storage. Allocated alongside weights%w
    ! when boole_delta_f is on. tau_c is the local collision time at the
    ! starting position; t_reg_on the switch-on time of the damping;
    ! nu_r the damping rate; trace_time the per-particle confinement time.
    real(dp), allocatable, public :: tau_c(:,:), t_reg_on(:,:)
    real(dp), allocatable, public :: nu_r(:,:),  trace_time(:,:)

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
    & boole_eliminate_particles_outside_flux, flux_threshold_for_elimination, boole_delta_f, &
    & profile_dir, equil_mapping_file, boole_constant_delta_B_psi, delta_B_psi_const, &
    & pert_m_fourier, pert_n_fourier, delta_B_psi_file, species_for_delta_f, &
    & nu_r_frac, n_collision_times_trace, m_collision_times_reg_on, coulomb_log, &
    & boole_compute_n_modes_dft, s_outer_cut, boole_skip_phase_for_test, &
    & s_inner_sample, s_outer_sample

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
! Allocates the per-particle regularisation arrays. Initialised to 0.
! ====================================================================
subroutine allocate_delta_f_per_particle(n_particles, n_species)

    integer, intent(in) :: n_particles, n_species

    if (.not. allocated(tau_c))      allocate(tau_c(n_particles, n_species))
    if (.not. allocated(t_reg_on))   allocate(t_reg_on(n_particles, n_species))
    if (.not. allocated(nu_r))       allocate(nu_r(n_particles, n_species))
    if (.not. allocated(trace_time)) allocate(trace_time(n_particles, n_species))

    tau_c      = 0.0_dp
    t_reg_on   = 0.0_dp
    nu_r       = 0.0_dp
    trace_time = 0.0_dp

end subroutine allocate_delta_f_per_particle

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
    real(dp)                                          :: vpar, vperp, t_tot, trace_time_n
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

    if (in%boole_delta_f) call allocate_delta_f_per_particle(n_particles, in%n_species)

    allocate(local_tetr_moments(moment_specs%n_moments,ntetr))
    kpart = 0
    iantithetic = 1
    if (in%boole_antithetic_variate) iantithetic = 2

    t_tot = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(counter, kpart, species, in, c, iantithetic, start, s, n_particles, trace_time) &
    !$OMP& REDUCTION(+:t_tot) &
    !$OMP& PRIVATE(p, l, n, i, x, vpar, vperp, t, ind_tetr, iface, local_tetr_moments, local_counter, particle_status, trace_time_n) &
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

            ! Initial trace-time bound. Overwritten on first orbit_timestep
            ! call from per-particle profile values when boole_delta_f.
            trace_time_n = start%t(species)

            i = 0

            do while (t%confined.lt.trace_time_n)
                i = i + 1

                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, t, x, vpar, vperp, ind_tetr, iface, species, iswmode_in=1)
                    t%step = t%step / start%v0(species)
                else
                    t%step = trace_time_n - t%confined
                endif

                call orbit_timestep_rmp_response_currents(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                                          local_tetr_moments, local_counter, species, trace_time_n)

                t%confined = t%confined + t%step - t%remain
                t_tot = t_tot + t%step - t%remain

                if (in%boole_delta_f .and. allocated(trace_time)) then
                    if (trace_time(n, species) > 0.0_dp) trace_time_n = trace_time(n, species)
                end if

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
        if (in%boole_delta_f) then
            ! Always start from w(0) = 0 (paper convention).
            weights%w(n, species) = (0.0_dp, 0.0_dp)
            weights%original(n, species) = (0.0_dp, 0.0_dp)
            ! Per-particle regularisation parameters from local profile values.
            call init_regularisation_for_particle(n, ind_tetr, x, vpar, vperp, species)
        end if
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

        ! Delta-f weight evolution (Albert 2016, Eq. 4) with linear
        ! regularisation switched on after t_reg_on. Updates weights%w
        ! in place using the exact ODE integrator over the dwell time.
        if (in%boole_delta_f .and. ind_tetr /= -1) then
            call update_delta_f_weight(n, ind_tetr, x, vpar, vperp, t_pass, t, species)
        endif

        call update_local_tetr_moments(local_tetr_moments, ind_tetr_save, n, optional_quantities, species)
        if ((grid_kind.eq.2).or.(grid_kind.eq.3)) call compute_radial_fluxes(ind_tetr_save, ind_tetr, x)

        ! Outer computational-domain cut: kill the particle on first entry
        ! into a tetra whose s exceeds s_outer_cut. Moments for the tetra just
        ! exited have already been accumulated above. Starting distribution
        ! and phase-space-volume normalisation are untouched.
        if (ind_tetr /= -1 .and. s_outer_cut < 1.0_dp) then
            if (eval_s_local(ind_tetr, x) > s_outer_cut) then
                ind_tetr = -1
                exit
            end if
        end if

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

    ! NOTE: the v_r(x_0) * f_0 multiplication that used to live in
    ! adapt_weights_delta_f is intentionally dropped. With the new scheme,
    ! the initial weight is the base PDF normalisation above, and the
    ! delta-f response accumulates along the orbit via update_delta_f_weight
    ! in orbit_timestep_rmp_response_currents.

end subroutine calc_particle_weights_and_jperp_rmp_response_currents

! ====================================================================
! Delta-f weight increment for a single tetrahedron crossing.
!
! Evaluates
!   wdot_s = -vpar * (dB_psi / B0) * (A1 + A2 * m v^2 / (2 T_alpha))
!   A1 = d(ln n)/dpsi_pol + (e/T) dPhi_0/dpsi_pol - (3/2) d(ln T_alpha)/dpsi_pol
!   A2 = d(ln T_alpha)/dpsi_pol
! at the new cell's entry point x, and accumulates w <- w + wdot_s * t_pass.
!
! The profiles (and their stored derivatives) are splined in
! s = psi_tor / psi_tor_edge. With q(s) = dpsi_tor/dpsi_pol, the chain
! rule gives
!
!   d/dpsi_pol = (dpsi_tor/dpsi_pol) * (ds/dpsi_tor) * d/ds
!              = q(s) / psi_tor_edge * d/ds.
!
! T_alpha is taken from Ti if the traced species is ionic, Te otherwise.
! ====================================================================
subroutine update_delta_f_weight(n, ind_tetr, x, vpar, vperp, t_pass, t, species)

    use gorilla_applets_types_mod, only: weights, time_t

    integer,      intent(in) :: n, ind_tetr, species
    real(dp),     intent(in) :: x(3), vpar, vperp, t_pass
    type(time_t), intent(in) :: t

    complex(dp) :: wdot_s
    real(dp)    :: t_current, nu_r_loc, decay
    complex(dp), parameter :: zero_c = (0.0_dp, 0.0_dp)

    call eval_wdot_s(ind_tetr, x, vpar, vperp, species, wdot_s)

    t_current = t%confined + t%step - t%remain

    if (t_current > t_reg_on(n, species)) then
        ! Regularised phase: dw/dt = wdot_s - nu_r * w. Exact integrator
        ! with wdot_s and nu_r constant over the dwell time t_pass.
        nu_r_loc = nu_r(n, species)
        if (nu_r_loc > 0.0_dp) then
            decay = exp(-nu_r_loc * t_pass)
            weights%w(n, species) = weights%w(n, species) * decay &
                + (wdot_s / cmplx(nu_r_loc, 0.0_dp, kind=dp)) * (1.0_dp - decay)
        else
            weights%w(n, species) = weights%w(n, species) + wdot_s * t_pass
        end if
    else
        ! Unregularised transient: pure dw/dt = wdot_s.
        weights%w(n, species) = weights%w(n, species) + wdot_s * t_pass
    end if

end subroutine update_delta_f_weight

! ====================================================================
! Evaluates the right-hand side wdot_s of the delta-f equation
! (Albert 2016, Eq. 4) at a given particle position and velocity.
! No state is mutated. Used both for the per-step time integration
! (update_delta_f_weight) and for the static-initial-weight option
! (set_initial_delta_f_weight).
! ====================================================================
subroutine eval_wdot_s(ind_tetr, x, vpar, vperp, species, wdot_s)

    use gorilla_applets_types_mod, only: in, start
    use tetra_physics_mod, only: tetra_physics, coord_system
    use tetra_grid_mod, only: tetra_grid, verts_sthetaphi
    use constants, only: ev2erg
    use profile_data_mod, only: eval_profiles, profile_values_t, eval_q, get_psi_tor_edge
    use perturbation_field_mod, only: eval_delta_B_psi

    integer,     intent(in)  :: ind_tetr, species
    real(dp),    intent(in)  :: x(3), vpar, vperp
    complex(dp), intent(out) :: wdot_s

    real(dp)    :: z_cell(3), s_loc, theta_loc, phi_loc
    real(dp)    :: B0_loc
    complex(dp) :: dB_psi
    real(dp)    :: v_sq, T_alpha_erg, A1, A2
    real(dp)    :: mass, charge, ds_dpsi_pol, q_loc, psi_tor_edge
    type(profile_values_t) :: pv

    wdot_s = (0.0_dp, 0.0_dp)
    z_cell = x - tetra_physics(ind_tetr)%x1

    s_loc = eval_s_local(ind_tetr, x)
    theta_loc = verts_sthetaphi(2, tetra_grid(ind_tetr)%ind_knot(1))

    if (coord_system .eq. 1) then
        phi_loc = x(2)
    else
        phi_loc = x(3)
    end if

    B0_loc = tetra_physics(ind_tetr)%bmod1 + sum(tetra_physics(ind_tetr)%gB * z_cell)
    if (B0_loc <= 0.0_dp) return

    call eval_delta_B_psi(s_loc, theta_loc, phi_loc, dB_psi)

    call eval_profiles(s_loc, pv)

    mass   = start%particle_mass(species)
    charge = start%particle_charge(species)
    v_sq   = vpar*vpar + vperp*vperp

    if (charge > 0.0_dp) then
        T_alpha_erg = pv%Ti * ev2erg
    else
        T_alpha_erg = pv%Te * ev2erg
    end if
    if (T_alpha_erg <= 0.0_dp) T_alpha_erg = in%energy_eV * ev2erg

    psi_tor_edge = get_psi_tor_edge()
    q_loc = eval_q(s_loc)
    if (abs(psi_tor_edge) > 0.0_dp .and. abs(q_loc) > 0.0_dp) then
        ds_dpsi_pol = q_loc / psi_tor_edge
    else
        ds_dpsi_pol = 0.0_dp
    end if

    if (charge > 0.0_dp) then
        A1 = (pv%dlnn_ds + (charge/T_alpha_erg) * pv%dPhi0_ds &
              - 1.5_dp * pv%dlnTi_ds) * ds_dpsi_pol
        A2 = pv%dlnTi_ds * ds_dpsi_pol
    else
        A1 = (pv%dlnn_ds + (charge/T_alpha_erg) * pv%dPhi0_ds &
              - 1.5_dp * pv%dlnTe_ds) * ds_dpsi_pol
        A2 = pv%dlnTe_ds * ds_dpsi_pol
    end if

    wdot_s = (-vpar * (A1 + A2 * mass * v_sq / (2.0_dp * T_alpha_erg))) &
             * (dB_psi / cmplx(B0_loc, 0.0_dp, kind=dp))

end subroutine eval_wdot_s

! ====================================================================
! Computes the local Maxwellian collision frequency for the chosen
! species at the starting position of particle n, then derives
! tau_c, t_reg_on, nu_r, and trace_time. Stored in module-level
! arrays for use by the per-step weight evolution and the trace-time
! cap.
! ====================================================================
subroutine init_regularisation_for_particle(n, ind_tetr, x, vpar, vperp, species)

    use gorilla_applets_types_mod, only: start
    use constants, only: pi, ev2erg
    use profile_data_mod, only: eval_profiles, profile_values_t

    integer,  intent(in) :: n, ind_tetr, species
    real(dp), intent(in) :: x(3), vpar, vperp

    real(dp) :: s_loc, n_loc, T_loc_erg, mass, charge, nu_c, tau_c_loc
    type(profile_values_t) :: pv

    s_loc = eval_s_local(ind_tetr, x)
    call eval_profiles(s_loc, pv)

    mass   = start%particle_mass(species)
    charge = start%particle_charge(species)
    n_loc  = max(pv%n_e, 1.0_dp)
    if (charge > 0.0_dp) then
        T_loc_erg = max(pv%Ti, 1.0_dp) * ev2erg
    else
        T_loc_erg = max(pv%Te, 1.0_dp) * ev2erg
    end if

    ! Spitzer-like collision frequency in CGS:
    !   nu_c = 4 pi n q^4 ln_Lambda / (sqrt(m) (2 T)^{3/2}).
    nu_c = 4.0_dp * pi * n_loc * charge**4 * coulomb_log &
         / (sqrt(mass) * (2.0_dp * T_loc_erg)**1.5_dp)
    if (nu_c <= 0.0_dp) nu_c = 1.0_dp
    tau_c_loc = 1.0_dp / nu_c

    tau_c(n, species)      = tau_c_loc
    nu_r(n, species)       = nu_r_frac * nu_c
    t_reg_on(n, species)   = real(m_collision_times_reg_on, dp) * tau_c_loc
    trace_time(n, species) = real(n_collision_times_trace,  dp) * tau_c_loc

end subroutine init_regularisation_for_particle

! ====================================================================
! Physics-normal flux label s at the particle position x inside tetra
! ind_tetr, clamped to [0,1] with 0 = magnetic axis and 1 = last closed
! flux surface. In flux coordinates (coord_system==2) GORILLA's x(1) is
! already s. In cylindrical coordinates (axisymmetric devices) we use
! normalised poloidal flux, inverted relative to the applet's exit-data
! convention (which reports 1 at the axis) so that the s_outer_cut
! namelist entry matches the physics-normal convention used in the
! delta-f design doc.
! ====================================================================
real(dp) function eval_s_local(ind_tetr, x) result(s_loc)

    use tetra_physics_mod, only: tetra_physics, coord_system
    use profile_data_mod, only: eval_s_from_psi_pol

    integer,  intent(in) :: ind_tetr
    real(dp), intent(in) :: x(3)

    real(dp) :: z_cell(3), psi_pol_loc

    z_cell = x - tetra_physics(ind_tetr)%x1

    if (coord_system == 2) then
        ! In symmetry flux coordinates GORILLA's x(1) is already s.
        s_loc = tetra_physics(ind_tetr)%x1(1) + z_cell(1)
    else
        ! Local poloidal flux (Aphi1 stores psi_pol on the reference vertex).
        psi_pol_loc = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi * z_cell)
        ! Convert to true s = psi_tor/psi_tor_edge using the spline built
        ! in profile_data_mod. This is independent of the mesh inner cut
        ! sfc_s_min, unlike the previous flux%poloidal_min/max version.
        s_loc = eval_s_from_psi_pol(psi_pol_loc)
    end if
    s_loc = max(0.0_dp, min(1.0_dp, s_loc))

end function eval_s_local

! ====================================================================
! Resamples the starting positions so they all land in the s-window
! [s_inner_sample, s_outer_sample].  For each particle, draws a fresh
! (R, phi, Z) (or (s, phi, theta)) from the same bounding box used by
! set_starting_positions, evaluates s, and accepts the draw if it falls
! in the window.  Up to max_tries redraws per particle; particles that
! cannot be placed are flagged lost.
! ====================================================================
subroutine bias_starting_positions_to_s_window()

    use gorilla_applets_types_mod, only: in, start, g
    use tetra_physics_mod, only: coord_system
    use find_tetra_mod, only: find_tetra
    use constants, only: pi

    integer, parameter :: max_tries = 1000
    integer  :: n, species, ind_tetr, iface, n_tries, n_replaced, n_failed
    real(dp) :: x(3), s_loc, u(3), x_amin, x_amax, x_cmin, x_cmax

    if (s_inner_sample <= 0.0_dp .and. s_outer_sample >= 1.0_dp) return

    print *, ''
    print *, 'Biasing start positions to s-window [', s_inner_sample, ',', s_outer_sample, ']'

    x_amin = g%amin; x_amax = g%amax
    x_cmin = g%cmin; x_cmax = g%cmax

    n_replaced = 0
    n_failed   = 0

    do species = 1, in%n_species
        do n = 1, in%num_particles
            if (start%lost(n, species)) cycle

            ! Quick check on the existing draw; keep it if already inside.
            x = start%x(:, n, species)
            call find_tetra(x, 0.0_dp, 0.0_dp, ind_tetr, iface)
            if (ind_tetr /= -1) then
                s_loc = eval_s_local(ind_tetr, x)
                if (s_loc >= s_inner_sample .and. s_loc <= s_outer_sample) cycle
            end if

            ! Redraw.
            do n_tries = 1, max_tries
                call random_number(u)
                if (coord_system .eq. 1) then
                    x(1) = x_amin + u(1) * (x_amax - x_amin)   ! R
                    x(2) = u(2) * 2.0_dp * pi                  ! phi
                    x(3) = x_cmin + u(3) * (x_cmax - x_cmin)   ! Z
                else
                    x(1) = x_amin + u(1) * (x_amax - x_amin)   ! s
                    x(2) = x_cmin + u(2) * (x_cmax - x_cmin)   ! theta
                    x(3) = u(3) * 2.0_dp * pi                  ! phi
                end if
                call find_tetra(x, 0.0_dp, 0.0_dp, ind_tetr, iface)
                if (ind_tetr == -1) cycle
                s_loc = eval_s_local(ind_tetr, x)
                if (s_loc >= s_inner_sample .and. s_loc <= s_outer_sample) then
                    start%x(:, n, species) = x
                    n_replaced = n_replaced + 1
                    exit
                end if
            end do
            if (n_tries > max_tries) then
                start%lost(n, species) = .true.
                n_failed = n_failed + 1
            end if
        end do
    end do

    print *, '  Replaced positions:', n_replaced
    print *, '  Could not place    :', n_failed
    print *, ''

end subroutine bias_starting_positions_to_s_window

end module utils_rmp_response_currents_mod
