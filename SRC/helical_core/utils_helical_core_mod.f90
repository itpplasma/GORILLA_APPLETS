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
               boole_monoenergetic, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, &
               boole_write_exit_data, boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, &
               boole_eliminate_particles_outside_flux, boole_delta_f, boole_unidirectional_beam
    integer :: i_integrator_type, seed_option, n_species, n_time_slices
    real(dp) :: flux_threshold_for_elimination, beam_density_fraction, reservoir_min_volume_fraction

    integer :: s_inp_unit

    NAMELIST /helical_core_nml/ time_step, energy_eV, n_particles, boole_squared_moments, boole_point_source, &
    & boole_collisions, boole_precalc_collisions, density, boole_refined_sqrt_g, boole_monoenergetic, &
    & boole_linear_density_simulation, boole_antithetic_variate, boole_linear_temperature_simulation, i_integrator_type, &
    & seed_option, boole_write_vertex_indices, boole_write_vertex_coordinates, boole_write_prism_volumes, &
    & boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_species, &
    & boole_eliminate_particles_outside_flux, flux_threshold_for_elimination, boole_delta_f, &
    & n_time_slices, boole_unidirectional_beam, beam_density_fraction, reservoir_min_volume_fraction

    !optional keys keep their input_t defaults when absent from the file
    n_time_slices = in%n_time_slices
    boole_unidirectional_beam = in%boole_unidirectional_beam
    beam_density_fraction = in%beam_density_fraction
    reservoir_min_volume_fraction = in%reservoir_min_volume_fraction

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
    in%n_time_slices = n_time_slices
    in%boole_unidirectional_beam = boole_unidirectional_beam
    in%beam_density_fraction = beam_density_fraction
    in%reservoir_min_volume_fraction = reservoir_min_volume_fraction

    ! Set defaults for anomalous transport-specific parameters (disabled)
    in%boole_calc_diffusion_coefficient = .false.
    in%i_scan_option = 0
    in%anomalous_diffusion_coefficient = 0.0_dp

    print *,'GORILLA_APPLETS: Loaded input data from helical_core.inp'

end subroutine read_helical_core_inp_into_type

! ====================================================================
subroutine parallelised_particle_pushing_helical_core(species, n_particles_in, boole_reuse_weights_in, i_slice_in)
!
! Parallelised particle pushing for helical core.
! Standard guiding-center orbit integration without anomalous displacement.
! Markers continue from start (update_start_type), so repeated calls push them in time slices;
! boole_reuse_weights_in keeps the weights of the previous slice instead of recomputing them, and
! i_slice_in gives every slice its own random number streams (the same streams in every slice would
! replay the same collision noise for each marker and heat the markers).
!
    use gorilla_applets_types_mod, only: counter, c, in, time_t, moment_specs, counter_t, particle_status_t, start, s
    use tetra_grid_mod, only: ntetr
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use utils_parallelised_particle_pushing_mod, only: print_progress, handle_lost_particles, add_marker_moments_to_thread_sums, &
        add_thread_moments_to_output, &
        add_local_counter_to_counter, initialise_loop_variables, carry_out_collisions, update_exit_data, update_start_type, &
        initialise_seed_for_random_numbers_for_each_thread
    use utils_orbit_timestep_mod, only: start_touched_tetr_tracking, stop_touched_tetr_tracking

    integer, intent(in)                               :: species
    integer, intent(in), optional                     :: n_particles_in
    logical, intent(in), optional                     :: boole_reuse_weights_in
    logical                                           :: boole_reuse_weights
    integer, intent(in), optional                     :: i_slice_in
    integer                                           :: i_slice
    integer                                           :: kpart, iantithetic, ind_tetr, iface, n_particles
    integer                                           :: p, l, n, i
    real(dp), dimension(3)                            :: x
    real(dp)                                          :: vpar, vperp, t_tot
    type(time_t)                                      :: t
    type(counter_t)                                   :: local_counter
    type(particle_status_t)                           :: particle_status
    complex(dp), dimension(:,:), allocatable          :: local_tetr_moments, thread_tetr_moments, thread_prism_moments_squared
    logical                                           :: thread_flag = .true.

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif
    boole_reuse_weights = .false.
    if (present(boole_reuse_weights_in)) boole_reuse_weights = boole_reuse_weights_in
    i_slice = 1
    if (present(i_slice_in)) i_slice = i_slice_in

    allocate(local_tetr_moments(moment_specs%n_moments,ntetr))
    kpart = 0
    iantithetic = 1
    if (in%boole_antithetic_variate) iantithetic = 2

    t_tot = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(counter, kpart, species, in, c, iantithetic, start, s, n_particles, boole_reuse_weights, i_slice, &
    !$OMP& moment_specs, ntetr) &
    !$OMP& REDUCTION(+:t_tot) &
    !$OMP& PRIVATE(p, l, n, i, x, vpar, vperp, t, ind_tetr, iface, local_tetr_moments, local_counter, particle_status, &
    !$OMP& thread_tetr_moments, thread_prism_moments_squared) &
    !$OMP& FIRSTPRIVATE(thread_flag)

    if (omp_get_thread_num().eq.0) print*, 'Number of threads: ', omp_get_num_threads()

    !moments are summed per thread and added to the output once per thread
    allocate(thread_tetr_moments(moment_specs%n_moments,ntetr))
    allocate(thread_prism_moments_squared(moment_specs%n_moments, &
                                          merge(ntetr/3, 0, moment_specs%boole_squared_moments)))
    thread_tetr_moments = (0.0_dp, 0.0_dp)
    thread_prism_moments_squared = (0.0_dp, 0.0_dp)
    local_tetr_moments = (0.0_dp, 0.0_dp)
    call start_touched_tetr_tracking(ntetr)

    !$OMP DO SCHEDULE(static)
    do p = 1, n_particles/iantithetic

        if ((.not.in%boole_precalc_collisions).and.thread_flag) then
            call initialise_seed_for_random_numbers_for_each_thread(omp_get_thread_num(), i_slice)
            thread_flag = .false.
        endif

        do l = 1, iantithetic
            n = (p-1)*iantithetic + l

            !$omp atomic update
            kpart = kpart + 1
            call print_progress(n_particles, kpart, n)

            ! Skip particles that were eliminated at initialization
            if (start%lost(n, species)) then
                call update_exit_data(.true., 0.0_dp, start%x(:,n,species), 0.0_dp, 0.0_dp, 0, n, species_in=species, ind_tetr=-1)
                cycle
            endif

            call initialise_loop_variables(l, n, local_counter, particle_status, t, local_tetr_moments, x, vpar, vperp, species)

            !zero time step only initialises the marker (tetrahedron and weight), so that the first collision
            !already uses the final particle weight
            t%step = 0.0_dp
            call orbit_timestep_helical_core(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                          local_tetr_moments, local_counter, species, start%t(species), boole_reuse_weights)
            t%step = start%t(species)

            i = 0

            do while (t%confined.lt.start%t(species))
                i = i + 1

                ! Apply collisions if enabled
                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, t, x, vpar, vperp, ind_tetr, iface, species,iswmode_in=1)
                    t%step = t%step / start%v0(species)
                else
                    t%step = start%t(species) - t%confined
                endif

                ! Perform guiding-center orbit integration (NO anomalous displacement)
                call orbit_timestep_helical_core(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                              local_tetr_moments, local_counter, species, start%t(species), &
                                              boole_reuse_weights)

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

        call add_marker_moments_to_thread_sums(local_tetr_moments, thread_tetr_moments, thread_prism_moments_squared)
    enddo
    !$OMP END DO

    !$omp critical
    call add_thread_moments_to_output(thread_tetr_moments, thread_prism_moments_squared, species)
    !$omp end critical
    call stop_touched_tetr_tracking
    !$OMP END PARALLEL

    print*, 'Total tracing time / number of particles: ', t_tot/n_particles, 's'

end subroutine parallelised_particle_pushing_helical_core

! ====================================================================
subroutine orbit_timestep_helical_core(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                    local_tetr_moments, local_counter, species, t_tot, boole_reuse_weights)
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
    use gorilla_applets_types_mod, only: counter_t, particle_status_t, start, in, time_t, g, weights
    use tetra_grid_settings_mod, only: grid_kind, sfc_s_min
    use utils_orbit_timestep_mod, only: initialize_constants_of_motion, compute_radial_fluxes, &
        identify_particles_entering_annulus, update_local_tetr_moments, calc_particle_weights_and_jperp

    real(dp), dimension(3), intent(inout)        :: x
    real(dp), intent(inout)                      :: vpar, vperp
    type(time_t), intent(inout)                  :: t
    type(particle_status_t), intent(inout)       :: particle_status
    integer, intent(inout)                       :: ind_tetr, iface
    integer, intent(in)                          :: n, species
    complex(dp), dimension(:,:), intent(inout)   :: local_tetr_moments
    type(counter_t), intent(inout)               :: local_counter
    real(dp), intent(in)                         :: t_tot
    logical, intent(in)                          :: boole_reuse_weights

    real(dp), dimension(3)                       :: z_save, x_new
    real(dp)                                     :: t_pass, perpinv
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
        if (.not.boole_reuse_weights) then
            call calc_particle_weights_and_jperp(n, z_save, vpar, vperp, ind_tetr, species)
            if (in%boole_delta_f) then
                call adapt_weights_delta_f(n, z_save, vpar, vperp, ind_tetr, species)
                weights%original(n, species) = weights%w(n, species)
            endif
        endif
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
                    !else
                        !print*, "particle pushing across the hole surrounding the magnetic axis was successful"
                    endif
                else
                    ! Particle left at the outer boundary - lost
                    exit
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

        ! Apply smooth fading of particle weights for t > t_tot*2/3 (only for delta-f method)
        if (in%boole_delta_f) call apply_weight_fading(n, species, t, t_tot)

        call update_local_tetr_moments(local_tetr_moments, ind_tetr_save, n, optional_quantities, species)
        if ((grid_kind.eq.2).or.(grid_kind.eq.3)) call compute_radial_fluxes(ind_tetr_save, ind_tetr, x)

        if (boole_t_finished) exit
    enddo

end subroutine orbit_timestep_helical_core

! ====================================================================
subroutine eliminate_particles_outside_flux_threshold
!
! Eliminates particles whose starting position has normalized poloidal flux
! greater than flux_threshold_for_elimination. Particles are marked as lost
! by setting start%lost to .true.
!
! The normalized poloidal flux is computed as:
!   s = (local_flux - flux_min) / (flux_max - flux_min)
! where s=0 at the magnetic axis and s=1 at the boundary.
!
    use gorilla_applets_types_mod, only: in, start
    use utils_data_pre_and_post_processing_mod, only: normalised_poloidal_flux
    use tetra_physics_mod, only: tetra_physics
    use find_tetra_mod, only: find_tetra

    integer  :: n, species, ind_tetr, iface, n_eliminated
    real(dp) :: vpar, vperp, local_flux, normalized_flux
    real(dp), dimension(3) :: x, z_local

    if (.not. in%boole_eliminate_particles_outside_flux) return

    n_eliminated = 0

    do species = 1, in%n_species
        do n = 1, in%num_particles
            if (start%lost(n, species)) cycle  ! Already marked as lost

            x = start%x(:, n, species)
            vpar = 0.0_dp
            vperp = 0.0_dp

            ! Find which tetrahedron this particle is in
            call find_tetra(x, vpar, vperp, ind_tetr, iface)

            if (ind_tetr == -1) then
                ! Particle is outside the grid - mark as lost
                start%lost(n, species) = .true.
                n_eliminated = n_eliminated + 1
                cycle
            endif

            ! Compute local poloidal flux
            z_local = x - tetra_physics(ind_tetr)%x1
            local_flux = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi * z_local)

            ! Compute normalized flux (0 at axis, 1 at boundary)
            normalized_flux = normalised_poloidal_flux(local_flux)

            ! Eliminate particle if flux exceeds threshold
            if (normalized_flux > in%flux_threshold_for_elimination) then
                start%lost(n, species) = .true.
                n_eliminated = n_eliminated + 1
            endif
        enddo
    enddo

    print *, 'Eliminated ', n_eliminated, ' particles with normalized flux > ', in%flux_threshold_for_elimination

end subroutine eliminate_particles_outside_flux_threshold

! ====================================================================
subroutine calc_v_r(z_save, vpar, vperp, ind_tetr, v_r)
!
! Computes the contravariant velocity component in the radial direction (v^r).
!
! The guiding-center velocity is computed from the polynomial pusher's
! first-order coefficient (dz/dtau = b + A*z), then converted to physical
! time and projected onto the gradient of s.
!
! Input:
!   z_save   - local position within tetrahedron (x - x1)
!   vpar     - parallel velocity
!   vperp    - perpendicular velocity
!   ind_tetr - tetrahedron index
!
! Output:
!   v_r      - contravariant velocity component in radial direction
!
    use gorilla_applets_types_mod, only: flux
    use tetra_physics_mod, only: tetra_physics, cm_over_e
    use constants, only: clight

    real(dp), dimension(3), intent(in) :: z_save
    real(dp), intent(in) :: vpar, vperp
    integer, intent(in) :: ind_tetr
    real(dp), intent(out) :: v_r

    real(dp) :: perpinv, k1, k3, bmod, phi_elec, vperp2, vpar2
    real(dp), dimension(4) :: z4, b, amat_in_z
    real(dp), dimension(4,4) :: amat
    real(dp), dimension(3) :: v_gc_tau, grad_r
    real(dp) :: dt_dtau

    ! Compute B-field magnitude at particle position
    bmod = tetra_physics(ind_tetr)%bmod1 + sum(tetra_physics(ind_tetr)%gb * z_save)

    ! Compute electric potential at particle position
    phi_elec = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi * z_save)

    ! Compute perpendicular invariant (negative by GORILLA convention)
    perpinv = -0.5_dp * vperp**2 / bmod

    ! Compute auxiliary quantities (matching pusher_tetra_poly convention)
    vperp2 = vperp**2
    vpar2 = vpar**2
    k1 = vperp2 + vpar2 + 2.0_dp * perpinv * tetra_physics(ind_tetr)%bmod1
    k3 = tetra_physics(ind_tetr)%Phi1 - phi_elec

    ! Build the 4-vector z = (z_save, vpar)
    z4(1:3) = z_save
    z4(4) = vpar

    ! Compute ODE coefficients b (inhomogeneous term)
    b(1:3) = (tetra_physics(ind_tetr)%curlh * k1 &
            + perpinv * tetra_physics(ind_tetr)%gBxh1) * cm_over_e &
            - clight * (2.0_dp * k3 * tetra_physics(ind_tetr)%curlh &
            + tetra_physics(ind_tetr)%gPhixh1)
    b(4) = perpinv * tetra_physics(ind_tetr)%gBxcurlA - clight / cm_over_e * tetra_physics(ind_tetr)%gPhixcurlA

    ! Compute ODE matrix A
    amat = 0.0_dp
    amat(1:3,1:3) = perpinv * cm_over_e * tetra_physics(ind_tetr)%alpmat &
                  - clight * tetra_physics(ind_tetr)%betmat
    amat(4,4) = perpinv * cm_over_e * tetra_physics(ind_tetr)%spalpmat &
              - clight * tetra_physics(ind_tetr)%spbetmat
    amat(1:3,4) = tetra_physics(ind_tetr)%curlA

    ! Compute A*z
    amat_in_z = matmul(amat, z4)

    ! Guiding-center velocity in tau: dz/dtau = b + A*z (first-order coefficient)
    v_gc_tau = b(1:3) + amat_in_z(1:3)

    ! Convert from tau to Hamiltonian time: v_gc = v_gc_tau / dt_dtau
    dt_dtau = tetra_physics(ind_tetr)%dt_dtau_const

    ! Compute gradient of r from gradient of A_phi normalized by poloidal flux at boundary
    grad_r = tetra_physics(ind_tetr)%gAphi / (flux%poloidal_max - flux%poloidal_min)

    ! Compute contravariant velocity component v^r = v_gc . grad_r
    v_r = sum(v_gc_tau * grad_r) / dt_dtau

end subroutine calc_v_r

! ====================================================================
subroutine adapt_weights_delta_f(n, z_save, vpar, vperp, ind_tetr, species)
!
! Sets particle weights for the delta-f method.
!
    use gorilla_applets_types_mod, only: start, in, g, weights
    use tetra_physics_mod, only: tetra_physics
    use volume_integrals_and_sqrt_g_mod, only: sqrt_g
    use supporting_functions_mod, only: bmod_func
    use constants, only: pi
    use marker_distribution_mod, only: linear_profile_slope

    integer, intent(in) :: n, ind_tetr, species
    real(dp), dimension(3), intent(in) :: z_save
    real(dp), intent(in) :: vpar, vperp

    real(dp) :: v_r, base_weight, r, z

    call calc_v_r(z_save, vpar, vperp, ind_tetr, v_r)

    !Densities follow the linear profile f = f_M*linear_profile_factor(r), so -v^r\partial f/\partial r turns into
    !-v^r*linear_profile_slope*f_M
    weights%w(n,species) = -weights%w(n,species) * v_r * linear_profile_slope

end subroutine adapt_weights_delta_f

! ====================================================================
subroutine apply_weight_fading(n, species, t, t_tot)
!
! Applies smooth fading of particle weights for t > t_tot*(1-alpha).
! weights%w = weights%original * 0.5*(1 + cos(pi*(t - t_tot*(1-alpha))/(t_tot*alpha)))
! This smoothly fades from original weight at t=t_tot*(1-alpha) to 0 at t=t_tot.
!
    use gorilla_applets_types_mod, only: weights, time_t
    use constants, only: pi

    integer, intent(in)       :: n, species
    type(time_t), intent(in)  :: t
    real(dp), intent(in)      :: t_tot

    real(dp) :: t_current, fade_factor, t_fade_start
    real(dp) :: alpha = 1.0_dp / 3.0_dp  ! Fraction of tracing time over which fading occurs

    ! Compute current time after this push: t%confined + (t%step - t%remain)
    t_current = t%confined + t%step - t%remain
    t_fade_start = t_tot * (1.0_dp - alpha)

    if (t_current > t_fade_start) then
        fade_factor = 0.5_dp * (1.0_dp + cos(pi * (t_current - t_fade_start) / (t_tot * alpha)))
        weights%w(n, species) = weights%original(n, species) * fade_factor
    endif

end subroutine apply_weight_fading

! ====================================================================
subroutine set_unidirectional_beam
!
! Turns the markers into a cold beam along the magnetic field: pitch +1 for all markers and a marker
! density of beam_density_fraction times the background density.
!
    use gorilla_applets_types_mod, only: in, start

    if (.not.in%boole_unidirectional_beam) return

    if (in%boole_antithetic_variate) then
        print*, 'Error: boole_unidirectional_beam requires boole_antithetic_variate = .false.'
        stop
    endif

    start%pitch = 1.0_dp
    start%density = in%density*in%beam_density_fraction

    print*, 'Unidirectional beam: pitch 1, marker density / background density = ', in%beam_density_fraction

end subroutine set_unidirectional_beam

! ====================================================================
subroutine write_collision_relaxation_diagnostics(file_unit, time, x, pitch, energy_eV, lost, u_b, temp_b)
!
! Writes one row of marker and background reservoir totals for the energy and momentum conserving
! collisions (species 1, first background ion species, cylindrical coordinates). Markers carry
! weights%w/n_particles physical particles, the reservoir of a tetrahedron dens_mat*tetr_volume
! particles; sums are in cgs, velocities are normalised by v0. Lost markers keep their state at the
! boundary (update_start_type), so the _lost columns are what left the domain.
! Parallel momentum P = sum m v_par is not an invariant in a torus (mirror force, varying b); the
! canonical toroidal angular momentum L = sum (m v_par h_phi + q/c A_phi) of the markers plus
! sum m_b u_par h_phi of the background is, up to the h_phi difference between a marker and the
! centroid of its tetrahedron. Background averages <.>_b are taken over the tetrahedra inside
! flux_threshold_for_elimination when markers outside it are eliminated, else over all.
!
! Columns: time [s], <v_par>_m/v0, <u_par>_b/v0, P_m, P_b, P_lost, E_m, E_b, E_lost, E_total,
!          L_m, L_b, L_lost, L_total, <T>_b [eV], T_m [eV] (kinetic energy in the frame moving
!          with <v_par>_m), min T_b [eV], max |u_par|_b/v0, number of confined markers
!
    use gorilla_applets_types_mod, only: in, c, start, weights
    use constants, only: ev2erg, clight, pi
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: tetra_grid, verts_rphiz
    use find_tetra_mod, only: find_tetra
    use utils_data_pre_and_post_processing_mod, only: normalised_poloidal_flux

    integer, intent(in)                  :: file_unit
    real(dp), intent(in)                 :: time
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:), intent(in)   :: pitch, energy_eV
    logical, dimension(:), intent(in)    :: lost
    real(dp), dimension(:), intent(in)   :: u_b, temp_b

    real(dp), dimension(size(pitch))     :: w_phys, w_lost, v, vpar, p_phi
    real(dp), dimension(size(c%tetr_volume)) :: n_b, h_phi_b
    logical, dimension(size(c%tetr_volume))  :: region
    real(dp), dimension(3) :: x_n, z, centroid
    real(dp), dimension(3,4) :: verts
    real(dp) :: m_t, q_t, m_b, n_m, p_m, e_m, n_bt, p_b, e_b, p_lost, e_lost, l_m, l_b, l_lost, &
                vpar_mean, temp_m, vpar_n, vperp_n
    integer :: n, i, k, ind_tetr, iface

    m_t = start%particle_mass(1)
    q_t = start%particle_charge(1)
    m_b = c%mass(1)

    w_phys = merge(0.0_dp, real(weights%w(:,1), dp)/in%n_particles, lost)
    w_lost = merge(real(weights%w(:,1), dp)/in%n_particles, 0.0_dp, lost)
    v = sqrt(2.0_dp*energy_eV*ev2erg/m_t)
    vpar = pitch*v

    !canonical toroidal angular momentum at the marker positions (zero weight where it cannot be located)
    p_phi = 0.0_dp
    do n = 1, size(pitch)
        if (real(weights%w(n,1), dp).eq.0.0_dp) cycle
        x_n = x(:,n)
        vpar_n = vpar(n)
        vperp_n = sqrt(max(v(n)**2 - vpar_n**2, 0.0_dp))
        call find_tetra(x_n, vpar_n, vperp_n, ind_tetr, iface)
        if (ind_tetr.eq.-1) cycle
        z = x_n - tetra_physics(ind_tetr)%x1
        p_phi(n) = m_t*vpar(n)*(tetra_physics(ind_tetr)%h2_1 + sum(tetra_physics(ind_tetr)%gh2*z)) + &
                   q_t/clight*(tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*z))
    enddo

    p_lost = m_t*sum(w_lost*vpar)
    e_lost = sum(w_lost*energy_eV)*ev2erg
    l_lost = sum(w_lost*p_phi)
    n_m = sum(w_phys)
    p_m = m_t*sum(w_phys*vpar)
    e_m = sum(w_phys*energy_eV)*ev2erg
    l_m = sum(w_phys*p_phi)
    vpar_mean = p_m/(m_t*n_m)
    temp_m = 2.0_dp/3.0_dp*(e_m/n_m - 0.5_dp*m_t*vpar_mean**2)/ev2erg

    do i = 1, size(c%tetr_volume)
        do k = 1, 4
            verts(:,k) = verts_rphiz(:,tetra_grid(i)%ind_knot(k))
        enddo
        verts(2,2:4) = verts(2,1) + modulo(verts(2,2:4) - verts(2,1) + pi, 2.0_dp*pi) - pi
        centroid = sum(verts, dim=2)/4.0_dp
        z = centroid - tetra_physics(i)%x1
        h_phi_b(i) = tetra_physics(i)%h2_1 + sum(tetra_physics(i)%gh2*z)
        region(i) = .true.
        if (in%boole_eliminate_particles_outside_flux) region(i) = normalised_poloidal_flux( &
            tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi*z)).le.in%flux_threshold_for_elimination
    enddo

    n_b = c%dens_mat(1,:)*c%tetr_volume
    p_b = m_b*sum(n_b*u_b)
    e_b = sum(n_b*(1.5_dp*temp_b*ev2erg + 0.5_dp*m_b*u_b**2))
    l_b = m_b*sum(n_b*u_b*h_phi_b)
    n_bt = sum(n_b, mask=region)

    write(file_unit,'(18ES17.9,I10)') time, vpar_mean/start%v0(1), &
        sum(n_b*u_b, mask=region)/n_bt/start%v0(1), p_m, p_b, p_lost, e_m, e_b, e_lost, e_m + e_b + e_lost, &
        l_m, l_b, l_lost, l_m + l_b + l_lost, sum(n_b*temp_b, mask=region)/n_bt, temp_m, &
        minval(temp_b), maxval(abs(u_b))/start%v0(1), count((.not.lost).and.(real(weights%w(:,1), dp).ne.0.0_dp))
    flush(file_unit)

end subroutine write_collision_relaxation_diagnostics

! ====================================================================
subroutine write_collision_relaxation_profiles
!
! Final state of the conserving collisions for radial profiles (cylindrical coordinates):
! background_reservoir.dat: per tetrahedron s, R, Z, physical volume [cm^3], density, u_par/v0, T [eV],
!                          reservoir cell (s, R, Z at the centroid, s = normalised poloidal flux)
! markers_final.dat:       per marker s, R, Z, v_par/v0, energy [eV], physical weight, lost flag
!
    use gorilla_applets_types_mod, only: in, c, start, weights
    use constants, only: pi, ev2erg
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: tetra_grid, verts_rphiz
    use find_tetra_mod, only: find_tetra
    use utils_data_pre_and_post_processing_mod, only: normalised_poloidal_flux

    real(dp), dimension(3,4) :: verts
    real(dp), dimension(3) :: centroid, x_n, z
    real(dp) :: s_n, v, vpar_n, vperp_n
    integer :: i, k, n, ind_tetr, iface, file_unit

    open(newunit=file_unit, file='background_reservoir.dat', status='replace')
    write(file_unit,'(A)') '#s R Z volume density u_par/v0 T cell'
    do i = 1, size(c%tetr_volume)
        do k = 1, 4
            verts(:,k) = verts_rphiz(:,tetra_grid(i)%ind_knot(k))
        enddo
        verts(2,2:4) = verts(2,1) + modulo(verts(2,2:4) - verts(2,1) + pi, 2.0_dp*pi) - pi
        centroid = sum(verts, dim=2)/4.0_dp
        z = centroid - tetra_physics(i)%x1
        s_n = normalised_poloidal_flux(tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi*z))
        write(file_unit,'(7ES17.9,I10)') s_n, centroid(1), centroid(3), c%tetr_volume(i), c%dens_mat(1,i), &
            c%vpar_mat(1,i)/start%v0(1), c%temp_mat(1,i), c%cell_of_tetr(i)
    enddo
    close(file_unit)

    open(newunit=file_unit, file='markers_final.dat', status='replace')
    write(file_unit,'(A)') '#s R Z v_par/v0 energy_eV weight lost'
    do n = 1, in%num_particles
        if (real(weights%w(n,1), dp).eq.0.0_dp) cycle
        x_n = start%x(:,n,1)
        v = sqrt(2.0_dp*start%energy(n,1)*ev2erg/start%particle_mass(1))
        vpar_n = start%pitch(n,1)*v
        vperp_n = sqrt(max(v**2 - vpar_n**2, 0.0_dp))
        s_n = -1.0_dp
        call find_tetra(x_n, vpar_n, vperp_n, ind_tetr, iface)
        if (ind_tetr.ne.-1) s_n = normalised_poloidal_flux(tetra_physics(ind_tetr)%Aphi1 + &
                                  sum(tetra_physics(ind_tetr)%gAphi*(x_n - tetra_physics(ind_tetr)%x1)))
        write(file_unit,'(6ES17.9,I3)') s_n, start%x(1,n,1), start%x(3,n,1), vpar_n/start%v0(1), start%energy(n,1), &
            real(weights%w(n,1), dp)/in%n_particles, merge(1, 0, start%lost(n,1))
    enddo
    close(file_unit)

end subroutine write_collision_relaxation_profiles

end module utils_helical_core_mod
