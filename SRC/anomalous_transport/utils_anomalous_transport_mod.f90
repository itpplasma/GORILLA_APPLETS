module utils_anomalous_transport_mod
!
! Module for anomalous transport calculations.
!
! Contains routines for:
!   - Reading anomalous transport input parameters
!   - Parallelised particle pushing (simplified, without self-consistent EF)
!   - Orbit timestep integration (simplified)
!
! Note: Displacement routines are in anomalous_transport_displacement_mod
!

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine read_anomalous_transport_inp_into_type

    use gorilla_applets_types_mod, only: in

    real(dp) :: time_step, energy_eV, n_particles, density, anomalous_diffusion_coefficient
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_monoenergetic, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, &
               boole_write_exit_data, boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, &
               boole_calc_diffusion_coefficient
    integer :: i_integrator_type, seed_option, n_species, i_scan_option

    integer :: s_inp_unit

    NAMELIST /anomalous_transport_nml/ time_step, energy_eV, n_particles, boole_squared_moments, boole_point_source, &
    & boole_collisions, boole_precalc_collisions, density, boole_refined_sqrt_g, boole_monoenergetic, &
    & boole_linear_density_simulation, boole_antithetic_variate, boole_linear_temperature_simulation, i_integrator_type, &
    & seed_option, boole_write_vertex_indices, boole_write_vertex_coordinates, boole_write_prism_volumes, &
    & boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_species, boole_calc_diffusion_coefficient, &
    & i_scan_option, anomalous_diffusion_coefficient

    open(newunit = s_inp_unit, file='anomalous_transport.inp', status='unknown')
    read(s_inp_unit,nml=anomalous_transport_nml)
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
    in%boole_calc_diffusion_coefficient = boole_calc_diffusion_coefficient
    in%i_scan_option = i_scan_option
    in%anomalous_diffusion_coefficient = anomalous_diffusion_coefficient

    print *,'GORILLA_APPLETS: Loaded input data from anomalous_transport.inp'

end subroutine read_anomalous_transport_inp_into_type

! ====================================================================
subroutine parallelised_particle_pushing_anomalous_transport(species, n_particles_in)
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
    use anomalous_transport_displacement_mod, only: anomalous_transport_displacement

    integer, intent(in)                               :: species
    integer, intent(in), optional                     :: n_particles_in
    integer                                           :: kpart, iantithetic, ind_tetr, iface, n_particles
    integer                                           :: p, l, n, i, k
    real(dp), dimension(3)                            :: x
    real(dp)                                          :: vpar, vperp, t_tot, t_step_s
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
    !$OMP& PRIVATE(p, l, n, i, k, x, vpar, vperp, t, ind_tetr, iface, local_tetr_moments, local_counter, particle_status, t_step_s) &
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
            if (in%boole_calc_diffusion_coefficient) then
                k = 1
                t_step_s = start%t(species) / s%k
            endif

            do while (t%confined.lt.start%t(species))
                i = i + 1

                ! Apply collisions if enabled
                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, t, x, vpar, vperp, ind_tetr, iface, species)
                    t%step = t%step / start%v0(species)
                else
                    t%step = t%step_anomalous_transport
                endif

                if (in%boole_calc_diffusion_coefficient) then
                    t_step_s = start%t(species)/s%k 
                    k = int(t%confined/t_step_s) + 1
                    t_step_s = t_step_s*k - t%confined
                endif

                ! Perform guiding-center orbit integration
                call orbit_timestep_anomalous_transport(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                              local_tetr_moments, local_counter, species, t_step_s, k)

                t%confined = t%confined + t%step - t%remain
                t_tot = t_tot + t%step - t%remain

                if (ind_tetr.ne.-1) call anomalous_transport_displacement(x, ind_tetr, iface, t%step, vpar, vperp)
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

end subroutine parallelised_particle_pushing_anomalous_transport

! ====================================================================
subroutine orbit_timestep_anomalous_transport(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                    local_tetr_moments, local_counter, species, t_step_s, k)
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
    use gorilla_applets_types_mod, only: counter_t, particle_status_t, start, in, time_t, g, s, flux, weights
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
    real(dp), intent(inout)                      :: t_step_s
    integer, intent(inout)                       :: k

    real(dp), dimension(3)                       :: z_save, x_new, z_local
    real(dp)                                     :: t_pass, perpinv, t_pusher, local_poloidal_flux, s_local
    real(dp)                                     :: dist_to_axis, rand_frac
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
                    !else
                        !print*, "particle pushing across the hole surrounding the magnetic axis was successful"
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

        t_pusher = t%remain
        if (in%boole_calc_diffusion_coefficient) t_pusher = min(t%remain, t_step_s)

        select case(ipusher)
            case(1)
                call pusher_tetra_rk(ind_tetr, iface, x, vpar, z_save, t%remain, t_pass, boole_t_finished, iper_phi)
            case(2)
                call pusher_tetra_poly(poly_order, ind_tetr, iface, x, vpar, z_save, t_pusher, &
                                       t_pass, boole_t_finished, iper_phi, optional_quantities)
        end select

        vperp = vperp_func(z_save, perpinv, ind_tetr_save)

        ! Record displacement statistics for diffusion coefficient calculation
        if (in%boole_calc_diffusion_coefficient) then
            if (boole_t_finished .and. (t%remain >= t_step_s-1.0d-20)) then
                if (t%remain > t_step_s) boole_t_finished = .false.
                t_step_s = start%t(species)/s%k + t_pass
                local_counter%tetr_pushings = local_counter%tetr_pushings - 1

                ! Compute normalized radial coordinate s from poloidal flux at current position
                ! Use current tetrahedron (ind_tetr) and recompute local coordinates
                z_local = x - tetra_physics(ind_tetr)%x1
                local_poloidal_flux = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*z_local)
                s_local = (local_poloidal_flux - flux%poloidal_min) / (flux%poloidal_max - flux%poloidal_min)

                !$omp critical
                s%delta_s(k) = s%delta_s(k) + (s_local - s%s0) * weights%w(n,species)
                s%delta_s_squared(k) = s%delta_s_squared(k) + (s_local - s%s0)**2 * weights%w(n,species)
                s%check(k) = s%check(k) + 1
                !$omp end critical

                k = k + 1
            endif
        endif

        t%remain = t%remain - t_pass
        t_step_s = t_step_s - t_pass

        call update_local_tetr_moments(local_tetr_moments, ind_tetr_save, n, optional_quantities, species)
        if ((grid_kind.eq.2).or.(grid_kind.eq.3)) call compute_radial_fluxes(ind_tetr_save, ind_tetr, x)

        if (boole_t_finished) exit
    enddo

end subroutine orbit_timestep_anomalous_transport

! ====================================================================
subroutine calc_diffusion_coefficient(filename_in, D_out)
!
! Calculates diffusion coefficient for anomalous transport validation.
! Initiates particles at a single point and tracks displacement statistics
! to verify that anomalous_transport_displacement produces the expected
! diffusion coefficient.
!
! In normalized flux coordinates (s): D_normalized = D_physical/r_minor^2
!
! Optional arguments:
!   filename_in - custom filename for output data (default: 'diffusion_coefficient_data.dat')
!   D_out       - returns the computed diffusion coefficient B
!
    use gorilla_applets_types_mod, only: in, s, start, weights
    use llsq_mod, only: llsq
    use utils_data_pre_and_post_processing_mod, only: initialize_exit_data, set_weights

    character(len=*), intent(in), optional :: filename_in
    real(dp), intent(out), optional :: D_out

    character(len=256) :: filename
    integer :: i, n_particles, file_id
    real(dp) :: A, B, offset, r_minor, D_physical, D_normalized
    real(dp), dimension(:), allocatable :: data_for_diffusion_coefficient
    real(dp), dimension(:,:,:), allocatable :: rand_matrix

    ! Set default filename or use provided one
    if (present(filename_in)) then
        filename = trim(filename_in)
    else
        filename = 'diffusion_coefficient_data.dat'
    endif

    ! Set number of test particles for diffusion coefficient calculation
    n_particles = 2000
    s%n_particles = n_particles

    ! Set up time sampling (1000 equally spaced intervals)
    s%k = 1000
    if (.not.allocated(s%delta_s)) allocate(s%delta_s(s%k))
    if (.not.allocated(s%delta_s_squared)) allocate(s%delta_s_squared(s%k))
    if (.not.allocated(s%time)) allocate(s%time(s%k))
    if (.not.allocated(s%check)) allocate(s%check(s%k))
    if (.not.allocated(data_for_diffusion_coefficient)) allocate(data_for_diffusion_coefficient(s%k))

    ! Initialize arrays
    s%delta_s = 0.0_dp
    s%delta_s_squared = 0.0_dp
    s%check = 0

    ! Reinitialize particles at a single point in cylindrical coordinates
    allocate(rand_matrix(5, n_particles, in%n_species))
    call RANDOM_NUMBER(rand_matrix)

    call allocate_start_type(n_particles)  ! Use local allocate_start_type with explicit particle count
    call allocate_weights(n_particles)
    call set_starting_positions_point_source(rand_matrix, (/1/))  ! This also sets s%s0
    call set_weights
    call set_rest_of_start_type(rand_matrix)  ! This also sets start%t(1)

    ! Initialize time array based on actual tracing time
    s%time = [(start%t(1)/s%k*i, i = 1, s%k)]

    deallocate(rand_matrix)

    ! Initialize exit data for tracking
    call initialize_exit_data(n_particles)

    print*, ''
    print*, '=============================================='
    print*, 'Diffusion Coefficient Calculation'
    print*, '=============================================='
    print*, 'Number of test particles: ', n_particles
    print*, 'Reference flux surface s0: ', s%s0
    print*, 'Number of time samples: ', s%k
    print*, 'Total tracing time: ', start%t(1), 's'
    print*, ''

    ! Call parallelised particle pushing (displacement tracking happens automatically)
    call parallelised_particle_pushing_anomalous_transport(species=1, n_particles_in=n_particles)

    ! Normalize by total weight
    do i = 1, s%k
        if (s%check(i) > 0) then
            s%delta_s(i) = s%delta_s(i) / sum(weights%w(:,1))
            s%delta_s_squared(i) = s%delta_s_squared(i) / sum(weights%w(:,1))
        endif
    enddo

    ! Perform least-squares fitting
    ! Discard first 20% of data to avoid transient effects
    i = ceiling(dble(s%k) * 0.2_dp)

    ! Fit convection coefficient A from <delta_s> vs time
    call llsq(int(s%k - i + 1, kind=8), s%time(i:s%k), s%delta_s(i:s%k), A, offset)

    ! Fit diffusion coefficient B from <(delta_s - A*t)^2> vs time
    data_for_diffusion_coefficient = s%delta_s_squared - 2*A*s%time*s%delta_s + A**2*s%time**2
    call llsq(int(s%k - i + 1, kind=8), s%time(i:s%k), data_for_diffusion_coefficient(i:s%k), B, offset)
    B = B / 2.0_dp

    ! Output results
    ! Expected value in normalized coordinates: D_normalized ~ D_physical / r_minor^2 s^-1
    r_minor = 50.0d0
    D_physical = in%anomalous_diffusion_coefficient
    D_normalized = D_physical/r_minor**2
    print*, '=============================================='
    print*, 'Diffusion Coefficient Results'
    print*, '=============================================='
    print*, 'Convection coefficient A: ', A, '1/s'
    print*, 'Diffusion coefficient B:  ', B, '1/s'
    print*, 'Expected value (normalized): ', D_normalized, '1/s'
    print*, 'Expected value (physical):   ', D_physical, 'cm^2/s'
    if (D_physical.gt.1.0d-20) print*, 'Relative error:           ', abs(B - D_normalized) / D_normalized * 100.0_dp, '%'
    print*, '=============================================='
    print*, ''

    ! Return diffusion coefficient if requested
    if (present(D_out)) D_out = B

    ! Write data to file (no comments for MATLAB compatibility)
    open(newunit=file_id, file=trim(filename))
    do i = 1, s%k
        write(file_id, '(4ES16.7)') s%time(i), s%delta_s(i), s%delta_s_squared(i), dble(s%check(i))
    enddo
    close(file_id)

    print*, 'Data written to: ', trim(filename)
    print*, ''

end subroutine calc_diffusion_coefficient

! ====================================================================
subroutine scan_anomalous_transport
!
! Unified scan subroutine for anomalous transport.
! Scans over different parameters based on i_scan_option:
!   1 = scan over eps_Phi (electric potential perturbation magnitude)
!   2 = scan over n3 (grid resolution in theta/Z direction)
!   3 = scan over n2 (grid resolution in phi direction)
!
! For grid scans (n2, n3), a fixed eps_Phi = 0.1 is used.
!
    use gorilla_applets_types_mod, only: in
    use tetra_grid_mod, only: nvert
    use tetra_grid_settings_mod, only: grid_kind, grid_size, n_extra_rings
    use tetra_physics_mod, only: phi_elec, make_tetra_physics
    use gorilla_settings_mod, only: coord_system
    use field_mod, only: ipert

    integer, parameter :: n_scan_points = 19
    real(dp), dimension(n_scan_points) :: scan_values, D_array
    real(dp), dimension(:), allocatable :: phi_elec_original
    real(dp) :: eps_Phi, D_coeff
    character(len=256) :: filename, summary_filename, scan_name, scan_var_name
    integer :: i_scan, file_id
    logical :: boole_grid_scan

    ! Check that anomalous_diffusion_coefficient is zero
    if (in%anomalous_diffusion_coefficient.gt.1.0d-10) then
        print*, ''
        print*, 'ERROR: anomalous_diffusion_coefficient must be set to zero'
        print*, '       when using i_scan_option > 0'
        print*, ''
        stop
    endif

    ! Setup scan parameters based on i_scan_option
    call setup_scan_parameters(in%i_scan_option, n_scan_points, scan_values, eps_Phi, &
                               summary_filename, scan_name, scan_var_name, boole_grid_scan)

    print*, ''
    print*, '=============================================='
    print*, 'Scanning Anomalous Transport over ', trim(scan_name)
    print*, '=============================================='
    print*, ''
    print*, 'Number of scan points: ', n_scan_points
    if (in%i_scan_option == 1) then
        print*, 'eps_Phi range: [', scan_values(1), ', ', scan_values(n_scan_points), ']'
    else
        print*, trim(scan_var_name), ' range: [', int(scan_values(1)), ', ', int(scan_values(n_scan_points)), ']'
        print*, 'Fixed eps_Phi: ', eps_Phi
    endif
    print*, 'Grid kind: ', grid_kind
    print*, ''

    ! Store original phi_elec (only needed for eps_Phi scan)
    if (.not.boole_grid_scan) then
        allocate(phi_elec_original(nvert))
        phi_elec_original = phi_elec
    endif

    ! Main scan loop
    do i_scan = 1, n_scan_points
        print*, '----------------------------------------------'
        print*, 'Scan point ', i_scan, ' of ', n_scan_points
        if (in%i_scan_option == 1) then
            print*, 'eps_Phi = ', scan_values(i_scan)
        else
            print*, trim(scan_var_name), ' = ', int(scan_values(i_scan))
        endif
        print*, '----------------------------------------------'

        ! Update scan variable and rebuild grid if needed
        call update_scan_variable(in%i_scan_option, scan_values(i_scan), eps_Phi, phi_elec_original)

        ! Apply electric potential perturbation
        call apply_phi_perturbation(eps_Phi)

        ! boole_keep_phi_elec=.true. prevents make_tetra_physics from overwriting phi_elec
        call make_tetra_physics(coord_system, ipert, boole_keep_phi_elec=.true.)

        ! Generate filename for this scan point
        select case(in%i_scan_option)
            case(1)
                write(filename, '(A,I3.3,A)') 'diffusion_coefficient_data_eps_Phi_', i_scan, '.dat'
            case(2)
                write(filename, '(A,I3.3,A)') 'diffusion_coefficient_data_n2_', i_scan, '.dat'
            case(3)
                write(filename, '(A,I3.3,A)') 'diffusion_coefficient_data_n3_', i_scan, '.dat'
        end select

        ! Calculate diffusion coefficient
        call calc_diffusion_coefficient(filename_in=trim(filename), D_out=D_coeff)

        D_array(i_scan) = D_coeff

        print*, 'Diffusion coefficient D = ', D_coeff
        print*, ''
    enddo

    ! Write summary file
    open(newunit=file_id, file=trim(summary_filename))
    do i_scan = 1, n_scan_points
        if (in%i_scan_option == 1) then
            write(file_id, '(2ES16.7)') scan_values(i_scan), D_array(i_scan)
        else
            write(file_id, '(I6, ES16.7)') int(scan_values(i_scan)), D_array(i_scan)
        endif
    enddo
    close(file_id)

    print*, '=============================================='
    print*, 'Scan Complete'
    print*, '=============================================='
    print*, 'Summary written to: ', trim(summary_filename)
    print*, ''
    print*, 'Results:'
    if (in%i_scan_option == 1) then
        print*, '  eps_Phi                  D_coeff'
    else
        print*, '  ', trim(scan_var_name), '       D_coeff'
    endif
    do i_scan = 1, n_scan_points
        if (in%i_scan_option == 1) then
            print*, '  ', scan_values(i_scan), '    ', D_array(i_scan)
        else
            print*, '  ', int(scan_values(i_scan)), '    ', D_array(i_scan)
        endif
    enddo
    print*, '=============================================='
    print*, ''

    if (allocated(phi_elec_original)) deallocate(phi_elec_original)

end subroutine scan_anomalous_transport

! ====================================================================
subroutine setup_scan_parameters(i_scan_option, n_scan_points, scan_values, eps_Phi, &
                                  summary_filename, scan_name, scan_var_name, boole_grid_scan)
!
! Sets up scan parameters based on i_scan_option.
!
    integer, intent(in) :: i_scan_option, n_scan_points
    real(dp), dimension(n_scan_points), intent(out) :: scan_values
    real(dp), intent(out) :: eps_Phi
    character(len=*), intent(out) :: summary_filename, scan_name, scan_var_name
    logical, intent(out) :: boole_grid_scan

    real(dp) :: val_min, val_max
    integer :: i

    select case(i_scan_option)
        case(1)
            ! Scan over eps_Phi
            scan_name = 'Phi Perturbation'
            scan_var_name = 'eps_Phi'
            summary_filename = 'diffusion_vs_eps_Phi.dat'
            boole_grid_scan = .false.
            val_min = 3.0d-2
            val_max = 3.0d-1
            eps_Phi = 0.0_dp  ! Will be set from scan_values in the loop
        case(2)
            ! Scan over n2
            scan_name = 'n2 (grid resolution)'
            scan_var_name = 'n2'
            summary_filename = 'diffusion_vs_n2.dat'
            boole_grid_scan = .true.
            val_min = 20.0_dp
            val_max = 200.0_dp
            eps_Phi = 0.1_dp
        case(3)
            ! Scan over n3
            scan_name = 'n3 (grid resolution)'
            scan_var_name = 'n3'
            summary_filename = 'diffusion_vs_n3.dat'
            boole_grid_scan = .true.
            val_min = 20.0_dp
            val_max = 200.0_dp
            eps_Phi = 0.1_dp
        case default
            print*, 'ERROR: Invalid i_scan_option = ', i_scan_option
            stop
    end select

    ! Generate linearly spaced scan values
    do i = 1, n_scan_points
        scan_values(i) = val_min + (val_max - val_min) * (i - 1.0_dp) / (n_scan_points - 1.0_dp)
    enddo

end subroutine setup_scan_parameters

! ====================================================================
subroutine update_scan_variable(i_scan_option, scan_value, eps_Phi, phi_elec_original)
!
! Updates the scan variable and rebuilds grid if needed.
! For eps_Phi scan: resets phi_elec to original values and updates eps_Phi
! For grid scans: sets n2 or n3 and rebuilds the entire grid
!
    use tetra_grid_mod, only: nvert, make_tetra_grid
    use tetra_grid_settings_mod, only: set_n2, set_n3
    use tetra_physics_mod, only: phi_elec, make_tetra_physics
    use tetra_physics_poly_precomp_mod, only: make_precomp_poly
    use gorilla_settings_mod, only: coord_system, boole_grid_for_find_tetra
    use field_mod, only: ipert
    use find_tetra_mod, only: grid_for_find_tetra
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use utils_data_pre_and_post_processing_mod, only: calc_poloidal_flux

    integer, intent(in) :: i_scan_option
    real(dp), intent(in) :: scan_value
    real(dp), intent(inout) :: eps_Phi
    real(dp), dimension(:), allocatable, intent(in) :: phi_elec_original

    select case(i_scan_option)
        case(1)
            ! Scan over eps_Phi: reset phi_elec and update eps_Phi
            eps_Phi = scan_value
            phi_elec = phi_elec_original

        case(2)
            ! Scan over n2: rebuild grid
            call set_n2(nint(scan_value))
            call rebuild_grid_and_physics()

        case(3)
            ! Scan over n3: rebuild grid
            call set_n3(nint(scan_value))
            call rebuild_grid_and_physics()
    end select

end subroutine update_scan_variable

! ====================================================================
subroutine rebuild_grid_and_physics()
!
! Rebuilds the tetrahedral grid and associated physics after changing grid parameters.
!
    use tetra_grid_mod, only: make_tetra_grid, verts_rphiz
    use tetra_physics_mod, only: make_tetra_physics
    use tetra_physics_poly_precomp_mod, only: make_precomp_poly
    use gorilla_settings_mod, only: coord_system, boole_grid_for_find_tetra
    use field_mod, only: ipert
    use find_tetra_mod, only: grid_for_find_tetra
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use utils_data_pre_and_post_processing_mod, only: calc_poloidal_flux, &
        set_moment_specifications, initialise_output, &
        calc_collision_coefficients_for_all_tetrahedra
    use gorilla_applets_types_mod, only: in

    print*, 'Rebuilding tetrahedral grid...'
    call make_tetra_grid()

    call make_tetra_physics(coord_system, ipert)
    print*, 'Physics calculation of mesh is finished'

    call make_precomp_poly()

    if (boole_grid_for_find_tetra) call grid_for_find_tetra

    ! Reinitialize output arrays and volume integrals with new grid dimensions
    call set_moment_specifications()
    call initialise_output()
    call calc_square_root_g()
    call calc_volume_integrals()
    call calc_poloidal_flux(verts_rphiz)

    ! Recalculate collision coefficients if collisions are enabled
    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra()

end subroutine rebuild_grid_and_physics

! ====================================================================
subroutine apply_phi_perturbation(eps_Phi)
!
! Applies non-axisymmetric perturbation to electrostatic potential.
! For field-aligned grids, vertices in extra rings near the O-point are excluded.
!
    use tetra_grid_mod, only: nvert
    use tetra_grid_settings_mod, only: grid_kind, grid_size, n_extra_rings
    use tetra_physics_mod, only: phi_elec

    real(dp), intent(in) :: eps_Phi

    real(dp), dimension(:), allocatable :: rnd_perturbation
    integer :: iv, iv_in_slice, n_verts_in_extra_rings

    allocate(rnd_perturbation(nvert))
    call random_number(rnd_perturbation)

    ! Compute number of vertices in extra rings (for field-aligned grids)
    if (grid_kind.eq.2 .or. grid_kind.eq.3) then
        n_verts_in_extra_rings = (1 + n_extra_rings) * grid_size(3)
    else
        n_verts_in_extra_rings = 0
    endif

    ! Apply perturbation
    do iv = 1, nvert
        if (grid_kind.eq.2 .or. grid_kind.eq.3) then
            iv_in_slice = modulo(iv-1, nvert/grid_size(2)) + 1
            if (iv_in_slice > n_verts_in_extra_rings) then
                phi_elec(iv) = phi_elec(iv) + eps_Phi * rnd_perturbation(iv)
            endif
        else
            phi_elec(iv) = phi_elec(iv) + eps_Phi * rnd_perturbation(iv)
        endif
    enddo

    deallocate(rnd_perturbation)

end subroutine apply_phi_perturbation

! ====================================================================
subroutine set_starting_positions_point_source(rand_matrix, species_in)
!
! Sets all particles at the same point in cylindrical coordinates (R, phi, Z).
! This point should be approximately at mid-radius in the computational domain.
! Also computes and sets s%s0 (the reference flux surface) from this position.
!
    use gorilla_applets_types_mod, only: in, start, g, s, flux
    use tetra_physics_mod, only: tetra_physics
    use find_tetra_mod, only: find_tetra

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix
    integer, dimension(:), intent(in), optional :: species_in

    integer, dimension(:), allocatable :: species
    integer :: i, ind_tetr, iface
    real(dp), dimension(3) :: x0, z_save
    real(dp) :: local_poloidal_flux, vpar_dummy, vperp_dummy

    if (present(species_in)) then
        allocate(species(size(species_in)))
        species = species_in
    else
        allocate(species(in%n_species))
        species = [(i, i=1, in%n_species)]
    endif

    ! Set all particles to the same point in cylindrical coordinates
    ! These values should be adjusted to match the actual device geometry
    ! For now, use approximate mid-radius values
    start%x(1,:,species) = 140.0_dp  ! R at mid-radius
    start%x(2,:,species) = 0.0_dp                       ! phi = 0
    start%x(3,:,species) = 0.0_dp  ! Z at midplane

    ! Compute s%s0 from the starting position
    x0 = start%x(:, 1, species(1))

    ! Find which tetrahedron contains this point
    ! Use representative velocity values (only geometric location matters)
    vpar_dummy = 1.0d6
    vperp_dummy = 1.0d6
    call find_tetra(x0, vpar_dummy, vperp_dummy, ind_tetr, iface)

    if (ind_tetr == -1) then
        print*, 'WARNING: Starting position is outside computational domain!'
        print*, 'Setting s%s0 = 0.5 as fallback'
        s%s0 = 0.5_dp
    else
        ! Compute local coordinates within tetrahedron
        z_save = x0 - tetra_physics(ind_tetr)%x1

        ! Compute poloidal flux at this position
        local_poloidal_flux = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi * z_save)

        ! Normalize to get s coordinate (no square root)
        s%s0 = (local_poloidal_flux - flux%poloidal_min) / (flux%poloidal_max - flux%poloidal_min)
    endif

    print*, 'Particles initialized at point: R =', start%x(1,1,species(1)), &
            ', phi =', start%x(2,1,species(1)), ', Z =', start%x(3,1,species(1))
    print*, 'Reference flux surface s0 =', s%s0

end subroutine set_starting_positions_point_source

! ====================================================================
subroutine set_rest_of_start_type(rand_matrix)
!
! Sets particle velocities, energies, and pitch angles.
! Copied from utils_data_pre_and_post_processing_mod to allow customization.
!
    use gorilla_applets_types_mod, only: in, start
    use tetra_physics_mod, only: cm_over_e, particle_charge, particle_mass
    use gorilla_applets_settings_mod, only: i_option
    use constants, only: echarge,ame,clight
    use constants, only: ev2erg

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix

    start%pitch(:,:) = 2*rand_matrix(4,:,:)-1 !pitch parameter
    start%energy = in%energy_eV
    if (.not. in%boole_monoenergetic) then
        start%energy = 5*in%energy_eV*rand_matrix(5,:,:) !boltzmann energy distribution
    endif

    if (in%boole_antithetic_variate) then
        start%x(:,1:in%num_particles:2,:) = start%x(:,2:in%num_particles:2,:)
        start%pitch(1:in%num_particles:2,:) = -start%pitch(2:in%num_particles:2,:)
        start%energy(1:in%num_particles:2,:) = start%energy(2:in%num_particles:2,:)
    endif

    start%particle_charge = particle_charge
    start%particle_mass = particle_mass
    start%cm_over_e = cm_over_e

    ! Set tracing time for diffusion coefficient calculation
    ! Use same value as in self-consistent EF example: 2.0 * 1.7e-4
    start%t(1) = 3.4d-4!1.0d-3

    start%v0 = sqrt(2.0_dp*in%energy_eV*ev2erg/start%particle_mass)

end subroutine set_rest_of_start_type

! ====================================================================
subroutine allocate_start_type(n_particles_in)
!
! Allocates arrays in the start type for particle properties.
! If n_particles_in is provided, uses that; otherwise uses in%num_particles.
!
    use gorilla_applets_types_mod, only: start, in

    integer, intent(in), optional :: n_particles_in
    integer :: n_particles

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

    ! Before allocating, deallocate if necessary
    call deallocate_start_type

    allocate(start%x(3, n_particles, in%n_species))
    allocate(start%pitch(n_particles, in%n_species))
    allocate(start%energy(n_particles, in%n_species))
    allocate(start%jperp(n_particles, in%n_species))
    allocate(start%lost(n_particles, in%n_species))
    allocate(start%particle_charge(in%n_species))
    allocate(start%particle_mass(in%n_species))
    allocate(start%cm_over_e(in%n_species))
    allocate(start%t(in%n_species))
    allocate(start%v0(in%n_species))

end subroutine allocate_start_type

! ====================================================================
subroutine allocate_weights(n_particles_in)
!
! Allocates the weights array.
! If n_particles_in is provided, uses that; otherwise uses in%num_particles.
!
    use gorilla_applets_types_mod, only: in, weights

    integer, intent(in), optional :: n_particles_in
    integer :: n_particles

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

    call deallocate_weights
    allocate(weights%w(n_particles, in%n_species))

end subroutine allocate_weights

! ====================================================================
subroutine deallocate_start_type
!
! Deallocates arrays in the start type.
!
    use gorilla_applets_types_mod, only: start

    if (allocated(start%x))               deallocate(start%x)
    if (allocated(start%pitch))           deallocate(start%pitch)
    if (allocated(start%energy))          deallocate(start%energy)
    if (allocated(start%jperp))           deallocate(start%jperp)
    if (allocated(start%lost))            deallocate(start%lost)
    if (allocated(start%particle_charge)) deallocate(start%particle_charge)
    if (allocated(start%particle_mass))   deallocate(start%particle_mass)
    if (allocated(start%cm_over_e))       deallocate(start%cm_over_e)
    if (allocated(start%t))               deallocate(start%t)
    if (allocated(start%v0))              deallocate(start%v0)

end subroutine deallocate_start_type

! ====================================================================
subroutine deallocate_weights
!
! Deallocates the weights array.
!
    use gorilla_applets_types_mod, only: weights

    if (allocated(weights%w)) deallocate(weights%w)

end subroutine deallocate_weights

end module utils_anomalous_transport_mod
