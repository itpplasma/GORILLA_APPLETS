module utils_scef_diffusion_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    !Ghost layer of the random walk (in%rw_ghost_layer_width > 0), per species, filled by calc_density_via_random_walk
    !and used by apply_rw_ghost_layer: shell densities of the source walkers (per shell volume, time_step and number of
    !walkers, before the exit-time correction), reservoir part of the shell densities per unit reservoir density
    !(absolute, from the inflow rate), and the vertex index of the wall
    real(dp), dimension(:,:), allocatable :: rw_source_shell, rw_ghost_shell
    integer, dimension(2) :: rw_ghost_vertex = 0

contains

subroutine calc_diffusion_coefficients(species)
    !Calibration run: pushes test particles of the given species and fits A, B
    !on every flux layer. Writes A_and_B_<species>.dat for later reuse.
    !With in%diffusion_scan_s0 > 0 only the flux layer closest to it is calibrated
    !and the result goes to A_and_B_scan_species<species>.dat instead.

    use gorilla_applets_types_mod, only: in, dc, start, s, g, exit_data
    use tetra_grid_settings_mod, only: grid_size
    use utils_data_pre_and_post_processing_mod, only: prepare_next_round_of_parallelised_particle_pushing, &
    calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared, initialize_exit_data
    use llsq_mod, only: llsq
    use tetra_physics_mod, only: particle_mass,tetra_physics
    use constants, only: ev2erg
    use russian_roulette_mod, only: prepare_russian_roulette
    use tetra_grid_settings_mod, only: sfc_s_min
    use utils_scef_particle_init_mod, only: allocate_start_type, set_particle_type_specifications, set_starting_positions, &
                                            set_rest_of_individual_particle_specifications, allocate_weights, &
                                            init_marker_distributions
    use utils_scef_particle_pushing_mod, only: parallelised_particle_pushing

    integer, intent(in) :: species

    integer :: ns, i, n_particles, file_id, n_ignored, ns_first, ns_last
    real(dp) :: extrapolation_factor, A, B, offset, tau_calibration, v0, v_max, weights_before_redistribution
    real(dp) :: standard_deviation, dummy
    integer :: pass, quality_unit, series_unit, i0, mid, n_lost, iq
    real(dp) :: cell_size, tau_new, tau_cap, mean_x, std_x, skewness, kurtosis, A_end, sigma_A_end, B_end, sigma_B_end, &
                offset_A, A_half1, A_half2, rms_residual, n_coll
    real(dp), dimension(5) :: q
    real(dp), dimension(5), parameter :: q_levels = (/0.01_dp, 0.05_dp, 0.5_dp, 0.95_dp, 0.99_dp/)
    real(dp), dimension(:), allocatable :: x_end, z_sorted
    real(dp), dimension(:,:,:), allocatable :: rand_matrix
    character(len=100) :: filename, ns_str, species_str
    real(dp), dimension(:), allocatable :: data_for_diffusion_coefficient, data_for_convection_coefficient, energy_start
    logical :: ignore_condition

    !Hardcoded per-species calibration duration. Refine for ions when needed.
    select case (species)
    case (1) !ions
        tau_calibration = 3.4d-2
    case (2) !electrons: 2 * tau_c_ei (NRL formula), preserves the previous behaviour
        tau_calibration = 3.4d-4
    case default
        print*, 'Error: calc_diffusion_coefficients called with unsupported species ', species
        stop
    end select
    if (in%diffusion_tau.gt.0.0_dp) tau_calibration = in%diffusion_tau

    s%n_particles = in%n_particles_diffusion
    call random_number(dummy)

    if (.not.allocated(rand_matrix)) allocate(rand_matrix(5,s%n_particles,1))
    if (.not.allocated(dc%A)) allocate(dc%A(grid_size(1)+1, 2)) !diffusion coefficients on every flux layer including boundaries, per species
    if (.not.allocated(dc%B)) allocate(dc%B(grid_size(1)+1, 2))
    if (.not.allocated(dc%grad_A)) allocate(dc%grad_A(grid_size(1), 2))
    if (.not.allocated(dc%grad_B)) allocate(dc%grad_B(grid_size(1), 2))
    s%k = 1000
    s%j = 100
    if (.not.allocated(s%delta_s)) allocate(s%delta_s(s%k))
    if (.not.allocated(s%delta_s_squared)) allocate(s%delta_s_squared(s%k))
    if (.not.allocated(s%time)) allocate(s%time(s%k))
    if (.not.allocated(s%f_v)) allocate(s%f_v(s%k,s%j))
    if (.not.allocated(s%check)) allocate(s%check(s%k))
    if (.not.allocated(data_for_diffusion_coefficient)) allocate(data_for_diffusion_coefficient(s%n_particles))
    if (.not.allocated(data_for_convection_coefficient)) allocate(data_for_convection_coefficient(s%n_particles))
    if (.not.allocated(energy_start)) allocate(energy_start(s%n_particles))
    if (.not.allocated(s%boole_large_distance)) allocate(s%boole_large_distance(s%n_particles))
    s%boole_large_distance = .false.
    if (allocated(s%boole_hit_axis)) deallocate(s%boole_hit_axis)
    allocate(s%boole_hit_axis(s%n_particles))
    if (.not.allocated(dc%frac_axis)) allocate(dc%frac_axis(grid_size(1)+1, 2), dc%frac_edge(grid_size(1)+1, 2))
    dc%frac_axis(:, species) = 1.0_dp   !boundary vertices are extrapolated, not calibrated
    dc%frac_edge(:, species) = 1.0_dp

    call allocate_start_type(s%n_particles)
    call allocate_weights(s%n_particles)
    call set_particle_type_specifications
    call init_marker_distributions(s%temperature)   !samples the calibration energies at s%temperature
    start%v0(species) = start%v0(species)*sqrt(s%temperature/in%energy_eV)
    call initialize_exit_data(s%n_particles)

    start%t(species) = tau_calibration !check afterwards if this was too little time

    s%time = [(start%t(species)/s%k*i,i = 1,s%k)]

    ns_first = 2
    ns_last = grid_size(1)
    if (in%diffusion_scan_s0.gt.0.0_dp) then
        ns_first = minloc(abs(g%s_vertices(2:grid_size(1))-in%diffusion_scan_s0), dim=1) + 1
        ns_last = ns_first
    endif

    !Quality of the calibration per flux surface (one line each): shape of the final displacement distribution,
    !the drift from the end positions with its statistical error, and the fit of <ds>(t)
    write(species_str, '(I0)') species
    open(newunit=quality_unit, file='calibration_quality_species'//trim(species_str)//'.dat', status='replace')
    write(quality_unit,'(a)') '# ns s0 tau[s] tau/local_collision_time N lost_fraction mean_ds std_ds sqrt(2*B_fit*tau) '// &
        'std_ds/cell_width skewness excess_kurtosis q01 q05 q50 q95 q99 (quantiles of (ds-mean)/std) A_fit A_end sigma_A_end '// &
        'A_fit_first_half A_fit_second_half rms_residual_of_<ds>(t)_fit B_fit B_end sigma_B_end n_passes '// &
        'fraction_mirrored_at_axis fraction_lost_at_edge (bias of A, B if either exceeds calibration_bias_tolerance)'
    close(quality_unit)

    do ns = ns_first,ns_last

        print*, 'ns = ', ns, '/', grid_size(1)


        !Initiate particles of `species` at the different flux surfaces leaving out the boundaries
        s%s0 = g%s_vertices(ns)
        cell_size = 0.5_dp*(g%s_vertices(ns+1) - g%s_vertices(ns-1))

        start%t(species) = tau_calibration
        if (in%diffusion_tau_collisions.gt.0.0_dp) start%t(species) = in%diffusion_tau_collisions*local_collision_time(s%s0, species)

    !Pass 2 only if the drift's statistical error exceeds in%diffusion_sigma_A_target: it scales like tau^(-1/2), so
    !tau is raised by (sigma_A/target)^2, capped so that the cloud stays within diffusion_max_spread_cells cell widths
    do pass = 1,2
        s%time = [(start%t(species)/s%k*i,i = 1,s%k)]
        print*, 'calibration tracing time = ', start%t(species), ' s'

        call RANDOM_NUMBER(rand_matrix)
        call set_starting_positions(rand_matrix,(/species/), s%s0)
        call set_rest_of_individual_particle_specifications(rand_matrix,boole_diffusion_coefficient_in = .true., &
                                                            species_in=(/species/),n_particles_in=s%n_particles)
        energy_start = start%energy(1:s%n_particles,species) !overwritten by the pushing


        call prepare_next_round_of_parallelised_particle_pushing(species)
        s%boole_hit_axis = .false.

        if ((ns.eq.ns_first).and.(pass.eq.1)) then
            call calc_collision_coefficients_for_all_tetrahedra(species)

            v0 = start%v0(species)
            v_max = v0*sqrt(start%epsilon_max)
            weights_before_redistribution = g%total_volume*in%density
            call prepare_russian_roulette(v0,v_max,weights_before_redistribution,10)
        endif

        ! open(23,file = 'energy_before_pushing.dat')
        !     do i = 1,s%n_particles
        !         write(23,*) start%energy(i,species)
        !     enddo
        ! close(23)
        call parallelised_particle_pushing(species = species,j = 1,boole_diffusion_coefficient = .true.,n_particles_in=s%n_particles)



        ! open(23,file = 'energy_after_pushing.dat')
        !     do i = 1,s%n_particles
        !         write(23,*) start%energy(i,species)
        !     enddo
        ! close(23)

        !Starting index (Throw away first 20 percent of values)
        i = ceiling(dble(s%k)*0.2_dp)


        call llsq ( int(s%k - i + 1, kind=8), s%time(i:s%k), s%delta_s(i:s%k), A, offset)
        offset_A = offset

        data_for_diffusion_coefficient = s%delta_s_squared-2*A*s%time*s%delta_s+A**2*s%time**2
        call llsq ( int(s%k - i + 1, kind=8), s%time(i:s%k), data_for_diffusion_coefficient, B, offset)
        B = B/2

        x_end = exit_data%x(1,1:s%n_particles,species)
        mean_x = sum(x_end)/s%n_particles
        std_x = sqrt(sum((x_end-mean_x)**2)/(s%n_particles-1))
        sigma_A_end = std_x/sqrt(dble(s%n_particles))/start%t(species)

        if ((pass.eq.2).or.(in%diffusion_sigma_A_target.le.0.0_dp)) exit
        if (sigma_A_end.le.in%diffusion_sigma_A_target) exit
        tau_new = start%t(species)*(sigma_A_end/in%diffusion_sigma_A_target)**2
        if ((in%diffusion_max_spread_cells.gt.0.0_dp).and.(std_x.gt.0.0_dp)) then
            tau_cap = start%t(species)*(in%diffusion_max_spread_cells*cell_size/std_x)**2
            if (tau_new.gt.tau_cap) then
                print*, 'calibration: sigma_A target not reachable within ', in%diffusion_max_spread_cells, &
                        ' cell widths of spread; tau capped (more particles needed)'
                tau_new = tau_cap
            endif
        endif
        if (tau_new.le.1.01_dp*start%t(species)) exit
        print*, 'calibration: sigma_A = ', sigma_A_end, ' > target ', in%diffusion_sigma_A_target, ', recalibrating'
        start%t(species) = tau_new
    enddo

        ! n_ignored = 0
        ! data_for_convection_coefficient = exit_data%x(1,:,species)-s%s0
        ! do i = 1,s%n_particles
        !     ignore_condition = abs(data_for_convection_coefficient(i)).gt.0.99_dp*min(s%s0-sfc_s_min,1.0_dp-s%s0)
        !     if (ignore_condition) then
        !         data_for_convection_coefficient(i) = 0.0_dp
        !         n_ignored = n_ignored + 1
        !     endif
        ! enddo

        ! A = sum(exit_data%x(1,:,species)-s%s0)/(s%n_particles*start%t(species))
        ! data_for_diffusion_coefficient = abs(exit_data%x(1,:,species)-(s%s0+A*start%t(species)))
        ! call quicksort(data_for_diffusion_coefficient, 1, s%n_particles)
        ! standard_deviation = data_for_diffusion_coefficient(int(s%n_particles*0.682689_dp))
        ! B = standard_deviation**2/(2*start%t(species))


        dc%A(ns, species) = A
        dc%B(ns, species) = B

        print*, 'A = ', A
        print*, 'B = ', B
    write(ns_str, '(I0)') ns
    write(species_str, '(I0)') species

    !Quality of this surface's calibration
    z_sorted = (x_end-mean_x)/std_x
    call quicksort(z_sorted, 1, s%n_particles)
    do iq = 1,5
        q(iq) = z_sorted(max(1,min(s%n_particles, nint(q_levels(iq)*s%n_particles))))
    enddo
    skewness = sum(z_sorted**3)/s%n_particles
    kurtosis = sum(z_sorted**4)/s%n_particles - 3.0_dp
    n_lost = count(exit_data%lost(1:s%n_particles,species).ne.0)
    dc%frac_axis(ns, species) = dble(count(s%boole_hit_axis(1:s%n_particles)))/s%n_particles
    dc%frac_edge(ns, species) = dble(n_lost)/s%n_particles
    if (dc%frac_axis(ns, species)+dc%frac_edge(ns, species).gt.in%calibration_bias_tolerance) &
        print*, 'calibration of surface ', ns, ' biased: fraction mirrored at the axis ', dc%frac_axis(ns, species), &
                ', lost at the edge ', dc%frac_edge(ns, species)
    A_end = (mean_x - s%s0)/start%t(species)
    B_end = sum((x_end-mean_x)**2)/s%n_particles/(2.0_dp*start%t(species))
    sigma_B_end = B_end*sqrt(sum(((x_end-mean_x)**2 - 2.0_dp*B_end*start%t(species))**2)/(s%n_particles-1)) &
                  /(2.0_dp*B_end*start%t(species))/sqrt(dble(s%n_particles))
    i0 = ceiling(dble(s%k)*0.2_dp)
    mid = (i0 + s%k)/2
    rms_residual = sqrt(sum((s%delta_s(i0:s%k) - (A*s%time(i0:s%k) + offset_A))**2)/(s%k - i0 + 1))
    call llsq ( int(mid - i0 + 1, kind=8), s%time(i0:mid), s%delta_s(i0:mid), A_half1, offset)
    call llsq ( int(s%k - mid, kind=8), s%time(mid+1:s%k), s%delta_s(mid+1:s%k), A_half2, offset)
    n_coll = start%t(species)/local_collision_time(s%s0, species)
    open(newunit=quality_unit, file='calibration_quality_species'//trim(species_str)//'.dat', position='append')
        write(quality_unit,'(i4,25es14.5,i3,2es14.5)') ns, s%s0, start%t(species), n_coll, dble(s%n_particles), &
            dble(n_lost)/s%n_particles, mean_x - s%s0, std_x, sqrt(2.0_dp*max(B,0.0_dp)*start%t(species)), std_x/cell_size, &
            skewness, kurtosis, q, A, A_end, sigma_A_end, A_half1, A_half2, rms_residual, B, B_end, sigma_B_end, pass, &
            dc%frac_axis(ns, species), dc%frac_edge(ns, species)
    close(quality_unit)

    !<ds>(t), <ds^2>(t) and the number of contributions per time sample, for checking the fits afterwards
    filename = 'calibration_series_species' // trim(species_str) // '_' // trim(ns_str) // '.dat'
    open(newunit=series_unit, file=filename)
        write(series_unit,'(a)') '# t[s] <ds> <ds^2> n_contributions (fit window: last 80 % of the samples)'
        do iq = 1,s%k
            write(series_unit,'(3es16.7,i9)') s%time(iq), s%delta_s(iq), s%delta_s_squared(iq), s%check(iq)
        enddo
    close(series_unit)

    filename = 'exit_s_values_species' // trim(species_str) // '_' // trim(ns_str) // '.dat'
    !columns: final s, start energy [eV], final energy [eV], lost (1) or not (0)
    open(newunit=file_id,file = filename)
        do i = 1,s%n_particles
            write(file_id,*) exit_data%x(1,i,species), energy_start(i), &
                0.5_dp*start%particle_mass(species)*(exit_data%vpar(i,species)**2+exit_data%vperp(i,species)**2)/ev2erg, &
                exit_data%lost(i,species)
        enddo
    close(file_id)


    ! filename = 'data' // trim(ns_str) // '.dat'
    ! open(23,file = filename)
    !     do i = 1,s%k
    !         write(23,*) s%time(i), s%delta_s(i), s%delta_s_squared(i)
    !     enddo
    ! close(23)

    ! open(23,file = 'f_v.dat')
    !     do i = 1,s%j
    !         write(23,*) s%f_v(:,i)
    !     enddo
    ! close(23)

    !print*, weights%w(:,species)

    enddo

    if (in%diffusion_scan_s0.gt.0.0_dp) then
        write(species_str, '(I0)') species
        filename = 'A_and_B_scan_species' // trim(species_str) // '.dat'
        open(newunit=file_id,file = filename)
            write(file_id,*) s%s0, dc%A(ns_first, species), dc%B(ns_first, species)
        close(file_id)
        return
    endif

    !Extrapolate values to the boundaries
    extrapolation_factor = (g%s_vertices(2)-g%s_vertices(1))/(g%s_vertices(3)-g%s_vertices(2))
    dc%A(1, species) =     dc%A(2, species) + (dc%A(2, species)-dc%A(3, species))*extrapolation_factor
    dc%B(1, species) = max(dc%B(2, species) + (dc%B(2, species)-dc%B(3, species))*extrapolation_factor, 0.0_dp)
    extrapolation_factor = (g%s_vertices(grid_size(1)+1)-g%s_vertices(grid_size(1)))/ &
                           (g%s_vertices(grid_size(1))-g%s_vertices(grid_size(1)-1))
    dc%A(grid_size(1)+1, species) =     dc%A(grid_size(1), species) &
                                      + (dc%A(grid_size(1), species)-dc%A(grid_size(1)-1, species))*extrapolation_factor
    dc%B(grid_size(1)+1, species) = max(dc%B(grid_size(1), species) &
                                      + (dc%B(grid_size(1), species)-dc%B(grid_size(1)-1, species))*extrapolation_factor, 0.0_dp)

    !Calculate gradients
    do ns = 1,grid_size(1)
        dc%grad_A(ns, species) = (dc%A(ns+1, species)-dc%A(ns, species))/(g%s_vertices(ns+1)-g%s_vertices(ns))
        dc%grad_B(ns, species) = (dc%B(ns+1, species)-dc%B(ns, species))/(g%s_vertices(ns+1)-g%s_vertices(ns))
    enddo

    write(species_str, '(I0)') species
    filename = 'A_and_B_species' // trim(species_str) // '.dat'
    open(newunit=file_id,file = filename)
    do i = 1,grid_size(1)+1
        !3rd: field -dPhi/ds [statV] of the calibration; 4th, 5th: fractions of test particles mirrored at the axis / lost at the edge
        write(file_id,*) dc%A(i, species), dc%B(i, species), vertex_electric_field(i), &
                         dc%frac_axis(i, species), dc%frac_edge(i, species)
    enddo
    close(file_id)

end subroutine calc_diffusion_coefficients

function local_collision_time(s0, species) result(tau)

    !Collision time [s] at the target n, T on flux surface s0 (NRL plasma formulary, cgs/eV):
    !electrons tau_e = 3.44e5 T_e^1.5/(n lambda_e), lambda_e = 24 - ln(sqrt(n)/T_e);
    !ions      tau_i = 2.09e7 T_i^1.5 sqrt(m_i/m_p)/(n Z^4 lambda_i), lambda_i = 23 - ln(Z^3 sqrt(2n)/T_i^1.5)

    use gorilla_applets_types_mod, only: start
    use constants, only: echarge, amp
    use scef_target_profiles_mod, only: target_density, target_temperature, rad_from_s

    real(dp), intent(in) :: s0
    integer, intent(in) :: species
    real(dp) :: tau, n, T, Z, lambda

    n = target_density(rad_from_s(s0))
    T = target_temperature(rad_from_s(s0), species)
    if (species.eq.1) then
        Z = start%particle_charge(1)/echarge
        n = n/Z   !target_density is the electron density
        lambda = 23.0_dp - log(Z**3*sqrt(2.0_dp*n)/T**1.5_dp)
        tau = 2.09d7*T**1.5_dp*sqrt(start%particle_mass(1)/amp)/(n*Z**4*lambda)
    else
        lambda = 24.0_dp - log(sqrt(n)/T)
        tau = 3.44d5*T**1.5_dp/(n*lambda)
    endif

end function local_collision_time

subroutine calc_density_via_random_walk(species)
    !Random-walk substitute for direct pushing of `species`. The field-free cached
    !A_and_B_<species>.dat is read in once per species; the convection coefficient is
    !then always derived from the current electric field (which reproduces the cached
    !one for E = 0), so the potential carries over between source steps.

    use gorilla_applets_types_mod, only: in, time_t, dc, ep, g, output, start, exit_data, weights, filenames, one_d
    use tetra_grid_settings_mod, only: grid_size, sfc_s_min
    use binsrc_mod, only: binsrc
    use tetra_grid_mod, only: ntetr
    use tetra_physics_mod, only: tetra_physics
    use find_tetra_mod, only: find_tetra
    use utils_scef_electric_potential_mod, only: fill_vector_parts_with_value
    use utils_orbit_timestep_mod, only: calc_particle_weights_and_jperp
    use scef_target_profiles_mod, only: target_temperature, rad_from_s

    integer, intent(in) :: species

    real(dp) :: delta_x, delta_t, xi, A, B, cell_size
    real(dp) :: s_wall, A_w, B_w, cell_w, dt_w, L_w, pos, ghost_time, vprime_w, inflow_per_nb
    real(dp), dimension(:), allocatable :: ghost_density
    integer :: k_wall, n_tries
    logical :: boole_ghost, boole_ghost_failed
    integer :: i, ns, k, num_steps_min, count_lost_particles, num_particles, particle_multiplication
    integer :: density_unit, n, ind_tetr, iface
    real(dp), dimension(3) :: x_particle, z_save
    real(dp) :: vpar_init, vperp_init
    real(dp), dimension(:), allocatable :: rw_density, rw_prism_densities
    real(dp), dimension(:), allocatable :: position, weight, exit_time
    real(dp), dimension(:), allocatable :: marker_shell   !density contribution per flux shell of one marker (mean over its walkers)
    type(time_t) :: t
    logical :: boole_lost, boole_use_fit_function   !set below (an initialiser would imply SAVE)

    particle_multiplication = max(int(1.0d4/in%num_particles),1)
    num_particles = in%num_particles*particle_multiplication

    allocate(rw_density(grid_size(1)), rw_prism_densities(ntetr/3))
    rw_density = 0.0_dp

    !>This should later be an input
    allocate(position(num_particles))
    allocate(weight(num_particles))
    allocate(exit_time(num_particles))

    !Physical weights: for ions apply the honest-tracing weight formula (sqrt(g)
    !Jacobian + Maxwell-Boltzmann) at each start position via
    !calc_particle_weights_and_jperp, exactly like the first orbit_timestep call
    !would do. Electrons copy the ion weights. Their positions already equal the ion start
    !positions (set_starting_positions); start%x(:,:,1) must not be copied here, because
    !honest ion tracing has overwritten it with the ion exit positions (update_start_type).
    if (species .eq. 1) then
        do n = 1, in%num_particles
            x_particle = start%x(:, n, 1)
            vpar_init = 0.0_dp
            vperp_init = 0.0_dp
            call find_tetra(x_particle, vpar_init, vperp_init, ind_tetr, iface)
            if (ind_tetr .eq. -1) then
                weights%w(n, 1) = 0.0_dp   !drop particles find_tetra can't locate, matches honest-tracing behaviour
                cycle
            endif
            z_save = x_particle - tetra_physics(ind_tetr)%x1
            call calc_particle_weights_and_jperp(n, z_save, vpar_init, vperp_init, ind_tetr, 1, &
                                                 boole_include_potential_in = .false.)
        enddo
    elseif ((species .eq. 2) .and. (.not. in%boole_static_ne)) then
        weights%w(:, 2)  = weights%w(:, 1)
    endif

    do i = 1,particle_multiplication
        position(i:num_particles:particle_multiplication) = start%x(1,:,species)
        weight(i:num_particles:particle_multiplication) = weights%w(:,species)
    enddo

    ! call RANDOM_NUMBER(position)
    ! call generate_distribution_one_minus_x2(position)
    ! allocate(dc%A(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    ! allocate(dc%B(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    ! allocate(dc%grad_A(grid_size(1)))
    ! allocate(dc%grad_B(grid_size(1)))
    ! allocate(g%s_vertices(grid_size(1)+1))
    ! g%s_vertices = verts_sthetaphi(1, [(grid_size(3)*(ns-1)+1,ns=1,grid_size(1)+1)])
    ! dc%A = 0.0_dp
    ! dc%B = 0.1_dp
    ! dc%grad_A = 0.0_dp
    ! dc%grad_B = 0.0_dp


    !With boole_measured_drift the walk uses the tabulated (measured) A and B, otherwise the fit of ln B
    boole_use_fit_function = .not.in%boole_measured_drift
    if (.not.dc%boole_cached_data_loaded(species)) then
        call get_diffusion_coefficient_data(species, boole_use_fit_function)
        dc%boole_cached_data_loaded(species) = .true.
    endif
    call calc_convection_coefficient_from_electric_field(species, boole_use_fit_function)

    t%step = in%time_step
    num_steps_min = 1000

    !Ghost layer: the outer layer s >= s_wall (snapped to the nearest vertex) is a reservoir at the target density.
    !Walkers are absorbed when they enter it (instead of at s = 1); the inflow from the reservoir is represented by
    !the ghost walkers below, whose amplitude apply_rw_ghost_layer fits after the exit-time correction.
    boole_ghost = in%rw_ghost_layer_width.gt.0.0_dp
    s_wall = 1.0_dp
    k_wall = grid_size(1)+1
    if (boole_ghost) then
        k_wall = minloc(abs(g%s_vertices - (1.0_dp - in%rw_ghost_layer_width)**2), 1)
        k_wall = min(max(k_wall, 3), grid_size(1))
        s_wall = g%s_vertices(k_wall)
        rw_ghost_vertex(species) = k_wall
    endif

    !Walkers i = (n-1)*particle_multiplication+1 ... n*particle_multiplication are copies of marker n. The sum
    !and sum of squares over markers of their per-shell contributions give the statistical error of the densities.
    allocate(marker_shell(grid_size(1)))
    if (.not.allocated(one_d%marker_sum)) allocate(one_d%marker_sum(grid_size(1),2), one_d%marker_sum_sq(grid_size(1),2))
    one_d%marker_sum(:,species) = 0.0_dp
    one_d%marker_sum_sq(:,species) = 0.0_dp

count_lost_particles = 0
    do i = 1,num_particles
        t%confined = 0.0_dp
        boole_lost = .false.
        if (boole_ghost) then
            if (position(i).ge.s_wall) then   !born in the reservoir: represented by the ghost walkers
                boole_lost = .true.
                count_lost_particles = count_lost_particles+1
            endif
        endif
        if (mod(i-1,particle_multiplication).eq.0) marker_shell = 0.0_dp
        do while ((t%confined.lt.t%step).and.(.not.boole_lost))

            call rw_coefficients(position(i), k, cell_size, A, B)

            delta_t = min(t%step/num_steps_min, cell_size**2/((abs(A)+abs(B))*4), t%step - t%confined) !control maximum possible jump

            rw_density(k) = rw_density(k) + delta_t*weight(i)
            marker_shell(k) = marker_shell(k) + delta_t*weight(i)
            t%confined = t%confined + delta_t

            call random_number(xi)
            xi = sqrt(12.0_dp)*(xi-0.5_dp)
            delta_x = sqrt(2*delta_t*max(B,0.0_dp))*xi + A*delta_t
            !print*, 'diffusive part over convective part = ', A,B,abs(sqrt(2*delta_t*B)/(A*delta_t)), &
            !sqrt(2*delta_t*B)*xi/(A*delta_t)

            position(i) = position(i) + delta_x

            if (position(i).lt.sfc_s_min) then
                position(i) = 2.0_dp*sfc_s_min - position(i)
            endif
            if ((position(i).gt.1.0_dp) .or. (boole_ghost.and.(position(i).ge.s_wall))) then
                boole_lost = .true.
                count_lost_particles = count_lost_particles+1
            endif
        enddo
        exit_time(i) = t%confined
        if (mod(i,particle_multiplication).eq.0) then
            marker_shell = marker_shell/particle_multiplication
            one_d%marker_sum(:,species) = one_d%marker_sum(:,species) + marker_shell
            one_d%marker_sum_sq(:,species) = one_d%marker_sum_sq(:,species) + marker_shell**2
        endif
    enddo

    !Scatter per-marker mean exit time back to exit_data%t_confined so the honest-tracing
    !accumulator in calc_average_charge_density_per_flux_layer works for the RW path too.
    !Each marker fans out to particle_multiplication trials; the mean over trials gives
    !the expected confinement time for that marker, matching the per-particle semantics
    !of the honest-tracing branch.
    do n = 1, in%num_particles
        exit_data%t_confined(n, species) = sum(exit_time((n-1)*particle_multiplication+1 : n*particle_multiplication)) &
                                          /particle_multiplication
    enddo

    if (count_lost_particles.lt.num_particles) print*, 'Warning: the tracing time (', t%step,'s) was so short that only ',&
                                                        count_lost_particles,'out of', num_particles, &
                                                        'particles of species ', species, ' left the computation domain'

    open(newunit=density_unit,file=filenames%density(species))
    do ns = 1, grid_size(1)
        write(density_unit,*) sfc_s_min + (1.0_dp - sfc_s_min) * (ns - 0.5_dp) / (grid_size(1)+1), rw_density(ns)
    enddo
    close(density_unit)

    do i = 1,grid_size(1)
        rw_density(i) = rw_density(i)/(g%s_shell_volumes(i)*t%step*num_particles)
        call fill_vector_parts_with_value(rw_prism_densities, g%prisms_per_flux_tube(i,:), rw_density(i))
    enddo

    if (boole_ghost) then
        if (.not.allocated(rw_source_shell)) allocate(rw_source_shell(grid_size(1),2), rw_ghost_shell(grid_size(1),2))
        rw_source_shell(:,species) = rw_density
        rw_ghost_shell(:,species) = 0.0_dp
        !Ghost walkers = reservoir particles that step into the inner region. The reservoir s >= s_wall has the density
        !n_b, i.e. n_b*V'(s_wall) particles per unit s. Start uniformly in [s_wall, s_wall + L_w] (L_w = largest possible
        !inward jump), one step dt_w with the coefficients at the wall, accepted if it lands inside: the accepted ones have
        !the entry depths of real incoming particles, and the acceptance f gives the inflow per step,
        !n_b*V'*L_w*f. Inflow rate Gamma = n_b*V'*L_w*f/dt_w; each of the N ghost walkers carries Gamma/N particles per
        !second, so the reservoir density of shell k is (Gamma/N)*sum(time in k)/V_k. They walk like the source walkers
        !until they are back in the layer or reach time_step.
        call rw_coefficients(s_wall - 1.0d-12, k, cell_w, A_w, B_w)
        boole_ghost_failed = B_w.le.0.0_dp
        if (.not.boole_ghost_failed) then
            allocate(ghost_density(grid_size(1)))
            ghost_density = 0.0_dp
            ghost_time = 0.0_dp
            n_tries = 0
            dt_w = min(t%step/num_steps_min, cell_w**2/((abs(A_w)+abs(B_w))*4))
            L_w = sqrt(6.0_dp*dt_w*B_w) + abs(A_w)*dt_w
            ghost_walkers: do i = 1,num_particles
                do
                    n_tries = n_tries+1
                    if (n_tries.gt.1000*num_particles) then
                        boole_ghost_failed = .true.
                        exit ghost_walkers
                    endif
                    call random_number(xi)
                    pos = s_wall + L_w*xi
                    call random_number(xi)
                    xi = sqrt(12.0_dp)*(xi-0.5_dp)
                    pos = pos + sqrt(2*dt_w*B_w)*xi + A_w*dt_w
                    if (pos.lt.s_wall) exit
                enddo
                t%confined = 0.0_dp
                boole_lost = .false.
                do while ((t%confined.lt.t%step).and.(.not.boole_lost))
                    call rw_coefficients(pos, k, cell_size, A, B)
                    delta_t = min(t%step/num_steps_min, cell_size**2/((abs(A)+abs(B))*4), t%step - t%confined)
                    ghost_density(k) = ghost_density(k) + delta_t
                    t%confined = t%confined + delta_t
                    call random_number(xi)
                    xi = sqrt(12.0_dp)*(xi-0.5_dp)
                    pos = pos + sqrt(2*delta_t*max(B,0.0_dp))*xi + A*delta_t
                    if (pos.lt.sfc_s_min) pos = 2.0_dp*sfc_s_min - pos
                    boole_lost = pos.ge.s_wall
                enddo
                ghost_time = ghost_time + t%confined
            enddo ghost_walkers
        endif
        if (boole_ghost_failed) then
            print*, 'Warning: ghost layer of species ', species, ': no inflow could be generated at s = ', s_wall, &
                    ' (B = ', B_w, ', A = ', A_w, '), the reservoir contributes nothing'
        else
            !V'(s_wall) from the real shell volumes on both sides of the wall
            vprime_w = 0.5_dp*(g%s_shell_volumes(k_wall-1)/(g%s_vertices(k_wall)-g%s_vertices(k_wall-1)) + &
                               g%s_shell_volumes(k_wall)/(g%s_vertices(k_wall+1)-g%s_vertices(k_wall)))
            inflow_per_nb = vprime_w*L_w*(real(num_particles,dp)/n_tries)/dt_w
            rw_ghost_shell(:,species) = inflow_per_nb*ghost_density/(num_particles*g%s_shell_volumes)
            print*, 'ghost layer of species ', species, ': wall at s = ', s_wall, ' (vertex ', k_wall, '), ', &
                    num_particles, ' ghost walkers, acceptance ', real(num_particles,dp)/n_tries, &
                    ', strip L = ', L_w, ', step ', dt_w, ' s, mean time inside ', ghost_time/num_particles, ' s'
        endif
    endif


    !This is a bit cumbersome having to use rw_prism_densities and a loop over all prisms instead of simply putting
    !output%prism_moments into fill_vector_parts_with_value, but that does not work because the latter does not accept
    !a complex argument. Think about a solution to this problem
    do i = 1,ntetr/3
        output%prism_moments(1,i,species) = complex(rw_prism_densities(i), 0.0_dp)
    enddo

    ep%mean_exit_time(species) = sum(exit_time)/num_particles
    print*, 'Total tracing time of species ', species, ' divided by number of particles is: ', &
    ep%mean_exit_time(species), 's'
    !print*, 'exit times are : ', exit_data%t_confined(:,species)

contains

    subroutine rw_coefficients(pos, k, cell_size, A, B)
        !Drift A and diffusion coefficient B of the walk at position pos (shell k, its width cell_size)
        real(dp), intent(in) :: pos
        integer, intent(out) :: k
        real(dp), intent(out) :: cell_size, A, B
        real(dp) :: B_c, s_c, s_lo, s_hi, c_thermo

        !binsrc finds k such that g%s_vertices(k-1) < pos < g%s_vertices(k)
        call binsrc(g%s_vertices,1,grid_size(1)+1,pos,k)
        k = k-1
        cell_size = g%s_vertices(k+1) - g%s_vertices(k)
        A = dc%A(k, species) + dc%grad_A(k, species)*(pos-g%s_vertices(k))
        B = dc%B(k, species) + dc%grad_B(k, species)*(pos-g%s_vertices(k))

        if (in%boole_measured_drift.and.(k.eq.1)) then
            !Axis cell: B -> 0 like s; the Ito drift B_2/s_2 replaces the measured one of the first
            !calibrated surface (whose cloud reaches the axis), the rest of the drift is taken from it
            B_c = dc%B(2, species)
            B = B_c*pos/g%s_vertices(2)
            A = B_c/g%s_vertices(2) + B*(dc%r_rest(2, species) + dc%qE_T(2, species))
        endif

        if (boole_use_fit_function) then
            !B from the fit of ln B; A holds q*E/T here, so the E_r drift uses the same B as the
            !diffusion. Below s_c, B -> 0 like s (B = (ds/dr)^2 D_r with finite D_r). Since the
            !volume element is constant in s, the walk needs the Ito drift dB/ds.
            s_c = g%s_vertices(4)
            if (pos.lt.s_c) then
                B_c = exp(sum(dc%polynomial_coefficients_for_B(:, species)*(/1.0_dp,s_c,s_c**2,0.0_dp/)))
                B = B_c*pos/s_c
                A = B*A + B_c/s_c
            else
                B = exp(sum(dc%polynomial_coefficients_for_B(:, species) &
                            *(/1.0_dp,pos,pos**2,0.0_dp/)))
                A = B*A + B*sum(dc%polynomial_coefficients_for_A(:, species) &
                                *(/1.0_dp,pos,0.0_dp/)) !dB/ds = B * d(ln B)/ds
            endif
        endif

        !Thermodiffusion: the walkers carry no energy, so the flux driven by the temperature gradient,
        !-c*B*n*d(ln T)/ds with c = D12/D11 - 3/2, enters as an explicit drift (central difference in s)
        !With thermodiffusion_table_file, c comes from the table at the current field (linear between vertices)
        c_thermo = in%thermodiffusion_coefficient
        if (dc%boole_thermo_table.and.(species.eq.2)) c_thermo = dc%c_thermo(k, species) &
            + (dc%c_thermo(k+1, species)-dc%c_thermo(k, species))*(pos-g%s_vertices(k))/cell_size
        if (c_thermo.ne.0.0_dp) then
            s_lo = max(pos - 1.0d-4, 0.0_dp)
            s_hi = min(pos + 1.0d-4, 1.0_dp)
            A = A - c_thermo*B &
                    *log(target_temperature(rad_from_s(s_hi), species)/target_temperature(rad_from_s(s_lo), species)) &
                    /(s_hi - s_lo)
        endif

    end subroutine rw_coefficients

end subroutine calc_density_via_random_walk

subroutine apply_rw_ghost_layer(species)

    !Ghost layer (in%rw_ghost_layer_width > 0), called after the exit-time correction. Inside the wall the density of a
    !random-walk species is n = c*n_source + n_b*n_reservoir (c = ep%exit_time_correction, n_reservoir = reservoir
    !part per unit reservoir density from calc_density_via_random_walk), beyond the wall the target. The reservoir
    !density n_b is the target at the wall (Dirichlet condition at s_wall); electrons use Z times the target, as in the
    !edge buffer.

    use gorilla_applets_types_mod, only: in, g, ep, output, start
    use tetra_grid_settings_mod, only: grid_size
    use tetra_grid_mod, only: ntetr
    use utils_scef_electric_potential_mod, only: fill_vector_parts_with_value
    use scef_target_profiles_mod, only: target_density, rad_from_s
    use constants, only: echarge

    integer, intent(in) :: species
    integer :: ns, k_wall, k_b, i
    real(dp) :: c, z_factor, n_b, value
    real(dp), dimension(:), allocatable :: n_target, prism_values

    k_wall = rw_ghost_vertex(species)
    if (k_wall.lt.3) return
    k_b = k_wall-1   !last shell inside the wall
    c = ep%exit_time_correction
    z_factor = 1.0_dp
    if (species.eq.in%n_species) z_factor = start%particle_charge(1)/echarge
    n_b = target_density(rad_from_s(g%s_vertices(k_wall)))*z_factor

    allocate(n_target(grid_size(1)), prism_values(ntetr/3))
    do ns = 1,grid_size(1)
        n_target(ns) = 0.5_dp*(target_density(rad_from_s(g%s_vertices(ns))) + &
                               target_density(rad_from_s(g%s_vertices(ns+1))))*z_factor
    enddo

    prism_values = 0.0_dp
    do ns = 1,grid_size(1)
        if (ns.lt.k_wall) then
            value = c*rw_source_shell(ns,species) + n_b*rw_ghost_shell(ns,species)
        else
            value = n_target(ns)
        endif
        call fill_vector_parts_with_value(prism_values, g%prisms_per_flux_tube(ns,:), value)
    enddo
    do i = 1,ntetr/3
        output%prism_moments(1,i,species) = complex(prism_values(i), 0.0_dp)
    enddo

    print*, 'ghost layer of species ', species, ': last inner shell: source part ', &
            c*rw_source_shell(k_b,species)/n_target(k_b), ', reservoir part ', &
            n_b*rw_ghost_shell(k_b,species)/n_target(k_b), ' of the target'
    print*, 'ghost layer: reservoir part / reservoir density per inner shell: ', rw_ghost_shell(1:k_b,species)

end subroutine apply_rw_ghost_layer

subroutine get_diffusion_coefficient_data(species, boole_use_fit_function)

    use gorilla_applets_types_mod, only: dc, g, in, start
    use tetra_grid_settings_mod, only: grid_size
    use utils_polyfit_mod, only: quadratic_fit
    use constants, only: ev2erg
    use scef_target_profiles_mod, only: target_temperature, rad_from_s

    integer, intent(in) :: species
    logical, intent(in) :: boole_use_fit_function

    character(len=100) :: filename, species_str
    integer :: id, ns, lower_bound, upper_bound, ios, n_v
    real(dp), dimension(:), allocatable :: x,y, qEcal_T
    logical :: success, boole_field_column
    real(dp) :: a0, a1, a2
    character(len=512) :: line


    if (.not.allocated(dc%A))                allocate(dc%A(grid_size(1)+1, 2))
    if (.not.allocated(dc%A_from_first_run)) allocate(dc%A_from_first_run(grid_size(1)+1, 2))
    if (.not.allocated(dc%B))                allocate(dc%B(grid_size(1)+1, 2))
    if (.not.allocated(dc%grad_A))           allocate(dc%grad_A(grid_size(1), 2))
    if (.not.allocated(dc%grad_B))           allocate(dc%grad_B(grid_size(1), 2))

    write(species_str, '(I0)') species
    filename = 'A_and_B_species' // trim(species_str) // '.dat'
    if (.not.allocated(dc%E_cal)) allocate(dc%E_cal(grid_size(1)+1, 2), dc%dB_meas(grid_size(1)+1, 2), &
                                            dc%r_rest(grid_size(1)+1, 2), dc%qE_T(grid_size(1)+1, 2))
    if (.not.allocated(dc%frac_axis)) allocate(dc%frac_axis(grid_size(1)+1, 2), dc%frac_edge(grid_size(1)+1, 2))
    dc%E_cal(:, species) = 0.0_dp
    boole_field_column = .true.
    open(newunit=id,file = filename,action='read')
    do ns = 1,grid_size(1)+1
        read(id,'(a)') line
        read(line,*,iostat=ios) dc%A(ns, species), dc%B(ns, species), dc%E_cal(ns, species), &
                                dc%frac_axis(ns, species), dc%frac_edge(ns, species)
        if (ios.ne.0) then   !files without the bias columns: bias unknown, taken as none
            dc%frac_axis(ns, species) = 0.0_dp
            dc%frac_edge(ns, species) = 0.0_dp
            read(line,*,iostat=ios) dc%A(ns, species), dc%B(ns, species), dc%E_cal(ns, species)
        endif
        if (ios.ne.0) then   !files of older versions: A and B only
            read(line,*) dc%A(ns, species), dc%B(ns, species)
            dc%E_cal(ns, species) = 0.0_dp
            boole_field_column = .false.
        endif
    enddo
    close(id)
    !the boundary vertices are always extrapolated
    dc%frac_axis(1, species) = 1.0_dp
    dc%frac_edge(grid_size(1)+1, species) = 1.0_dp
    if (in%boole_measured_drift.and.(.not.boole_field_column)) print*, 'Warning: ', trim(filename), &
        ' has no calibration-field column; the calibration is assumed to be field-free'

    if (in%boole_measured_drift) then
        if ((in%thermodiffusion_coefficient.ne.0.0_dp).or.(len_trim(in%thermodiffusion_table_file).gt.0)) print*, &
            'Warning: boole_measured_drift together with a thermodiffusion model: the measured drift already contains it'
        !B' of the tabulated B (central differences, one-sided at the ends) and the measured rest of the drift,
        !r = (A - B')/B - qE_cal/T, which sets the zero-flux shape n'/n = r + qE/T of the walk
        n_v = grid_size(1)+1
        allocate(qEcal_T(n_v))
        do ns = 1,n_v
            qEcal_T(ns) = start%particle_charge(species)*dc%E_cal(ns, species) &
                          /(target_temperature(rad_from_s(g%s_vertices(ns)), species)*ev2erg)
        enddo
        do ns = 1,n_v
            if (ns.eq.1) then
                dc%dB_meas(ns, species) = (dc%B(2, species)-dc%B(1, species))/(g%s_vertices(2)-g%s_vertices(1))
            elseif (ns.eq.n_v) then
                dc%dB_meas(ns, species) = (dc%B(n_v, species)-dc%B(n_v-1, species))/(g%s_vertices(n_v)-g%s_vertices(n_v-1))
            else
                dc%dB_meas(ns, species) = (dc%B(ns+1, species)-dc%B(ns-1, species))/(g%s_vertices(ns+1)-g%s_vertices(ns-1))
            endif
            if (dc%B(ns, species).gt.0.0_dp) then
                dc%r_rest(ns, species) = (dc%A(ns, species) - dc%dB_meas(ns, species))/dc%B(ns, species) &
                    - qEcal_T(ns)
            else
                dc%r_rest(ns, species) = 0.0_dp
            endif
        enddo
        if (dc%B(n_v, species).le.0.0_dp) dc%r_rest(n_v, species) = dc%r_rest(n_v-1, species)
        print*, 'measured drift: r = (A-dB/ds)/B - qE_cal/T at the vertices: ', dc%r_rest(:, species)
        if (in%drift_smoothing_width.gt.0.0_dp) call smooth_measured_drift(species, qEcal_T)
    endif

    dc%A_from_first_run(:, species) = dc%A(:, species)
    if (boole_use_fit_function) dc%A(:, species) = 0.0_dp

    if ((species.eq.2).and.(len_trim(in%thermodiffusion_table_file).gt.0)) call read_thermodiffusion_table()
    if (.not.allocated(dc%c_thermo)) allocate(dc%c_thermo(grid_size(1)+1, 2))
    dc%c_thermo(:, species) = 0.0_dp

    do ns = 1,grid_size(1)
        dc%grad_A(ns, species) = (dc%A(ns+1, species)-dc%A(ns, species))/(g%s_vertices(ns+1)-g%s_vertices(ns))
        dc%grad_B(ns, species) = (dc%B(ns+1, species)-dc%B(ns, species))/(g%s_vertices(ns+1)-g%s_vertices(ns))
    enddo



    lower_bound = 2
    upper_bound = grid_size(1) - 4

    !Fit ln B (not B): B varies by orders of magnitude with the temperature profile, and a fit of B
    !itself can turn negative towards the edge. Surfaces with B <= 0 (noise) are left out.
    x = pack(g%s_vertices(lower_bound:upper_bound), dc%B(lower_bound:upper_bound, species) > 0.0_dp)
    y = log(pack(dc%B(lower_bound:upper_bound, species), dc%B(lower_bound:upper_bound, species) > 0.0_dp))

    call quadratic_fit(size(x), x, y, a0, a1, a2, success) ! Least-squares fit of ln B ≈ a2*s^2 + a1*s + a0

    dc%polynomial_coefficients_for_B(:, species) = (/a0,a1,a2,0.0_dp/)   !of ln B
    dc%polynomial_coefficients_for_A(:, species) = (/a1,2*a2,0.0_dp/)    !of d(ln B)/ds

end subroutine get_diffusion_coefficient_data

subroutine smooth_measured_drift(species, qEcal_T)

    !Smooths the measured A and B directly and only then forms the drift rest r: A, B and dB/ds at every
    !vertex are replaced by a local linear least-squares fit with Gaussian weights of width
    !in%drift_smoothing_width in s (dB/ds = slope of the fit of B), and r = (A - dB/ds)/B - qE_cal/T is
    !computed from the smoothed values. Smoothing r itself amplifies the noise of A where B is small (edge
    !electrons). With this r the walk drift is A_smooth + B_smooth*(qE/T - qE_cal/T), i.e. no division by B
    !survives. Only unbiased surfaces enter the fit: calibrated (not the extrapolated boundary vertices),
    !B > 0, and at most in%calibration_bias_tolerance of the test particles mirrored at the axis or lost at
    !the edge. Biased vertices get the value extrapolated from their unbiased neighbours. Where the smoothed
    !B is not positive, B = 0, r = 0 and dB/ds is set to the smoothed A, so that the drift is still A_smooth.

    use gorilla_applets_types_mod, only: dc, g, in
    use tetra_grid_settings_mod, only: grid_size

    integer, intent(in) :: species
    real(dp), dimension(:), intent(in) :: qEcal_T
    integer :: n_v, ns
    logical, dimension(:), allocatable :: good
    real(dp), dimension(:), allocatable :: w, d, A_new, B_new, dB_new
    real(dp) :: S0, S1, S2, det

    n_v = grid_size(1)+1
    allocate(good(n_v), w(n_v), d(n_v), A_new(n_v), B_new(n_v), dB_new(n_v))
    good = (dc%frac_axis(:, species) + dc%frac_edge(:, species).le.in%calibration_bias_tolerance) &
           .and.(dc%B(:, species).gt.0.0_dp)
    good(1) = .false.
    good(n_v) = .false.
    if (count(good).lt.3) then
        print*, 'Warning: fewer than 3 unbiased calibrated surfaces, no smoothing of the measured drift'
        return
    endif

    do ns = 1,n_v
        d = g%s_vertices - g%s_vertices(ns)
        w = merge(exp(-0.5_dp*(d/in%drift_smoothing_width)**2), 0.0_dp, good)
        S0 = sum(w)
        S1 = sum(w*d)
        S2 = sum(w*d**2)
        det = S0*S2 - S1**2
        if (det.gt.1.0d-12*S0*S2) then
            A_new(ns)  = (S2*sum(w*dc%A(:, species)) - S1*sum(w*d*dc%A(:, species)))/det
            B_new(ns)  = (S2*sum(w*dc%B(:, species)) - S1*sum(w*d*dc%B(:, species)))/det
            dB_new(ns) = (S0*sum(w*d*dc%B(:, species)) - S1*sum(w*dc%B(:, species)))/det
        else   !only one surface with noticeable weight: its values, no slope
            A_new(ns)  = sum(w*dc%A(:, species))/S0
            B_new(ns)  = sum(w*dc%B(:, species))/S0
            dB_new(ns) = 0.0_dp
        endif
    enddo

    do ns = 1,n_v
        if (B_new(ns).gt.0.0_dp) then
            dc%B(ns, species) = B_new(ns)
            dc%dB_meas(ns, species) = dB_new(ns)
            dc%r_rest(ns, species) = (A_new(ns) - dB_new(ns))/B_new(ns) - qEcal_T(ns)
        else
            dc%B(ns, species) = 0.0_dp
            dc%dB_meas(ns, species) = A_new(ns)
            dc%r_rest(ns, species) = 0.0_dp
        endif
    enddo
    print*, 'drift smoothing of A and B (width ', in%drift_smoothing_width, ' in s): ', count(good), ' of ', n_v, &
            ' vertices unbiased; biased: ', pack([(ns, ns=1,n_v)], .not.good)
    print*, 'smoothed A at the vertices: ', A_new
    print*, 'smoothed B at the vertices: ', dc%B(:, species)
    print*, 'smoothed dB/ds at the vertices: ', dc%dB_meas(:, species)
    print*, 'drift r from the smoothed A and B at the vertices: ', dc%r_rest(:, species)

end subroutine smooth_measured_drift

subroutine calc_convection_coefficient_from_electric_field(species, boole_use_fit_function)

    use tetra_physics_mod, only: tetra_physics
    use constants, only: ev2erg
    use gorilla_applets_types_mod, only: dc, in, s, start, g
    use tetra_grid_settings_mod, only: grid_size
    use scef_target_profiles_mod, only: target_temperature, rad_from_s

    integer, intent(in) :: species
    logical, intent(in) :: boole_use_fit_function
    integer :: i, ns
    real(dp) :: electric_field

    do i = 1,grid_size(1)+1
        electric_field = vertex_electric_field(i)

        !q*E/T; with the fit function the random walk multiplies it by its fitted B, otherwise by the
        !tabulated B here
        dc%A(i, species) = start%particle_charge(species)*electric_field &
                           /(target_temperature(rad_from_s(g%s_vertices(i)), species)*ev2erg)
        if (in%boole_measured_drift) then
            !measured drift with the calibration field replaced by the current one
            dc%qE_T(i, species) = dc%A(i, species)
            dc%A(i, species) = dc%dB_meas(i, species) + dc%B(i, species)*(dc%r_rest(i, species) + dc%qE_T(i, species))
        elseif (.not.boole_use_fit_function) then
            dc%A(i, species) = dc%B(i, species)*dc%A(i, species) + dc%A_from_first_run(i, species)
        endif
        if (dc%boole_thermo_table.and.(species.eq.2)) &
            dc%c_thermo(i, species) = thermodiffusion_table_value(g%s_vertices(i), electric_field)
    enddo
    if (dc%boole_thermo_table.and.(species.eq.2)) then
        print*, 'thermodiffusion c from the table at the vertices: ', dc%c_thermo(:, species)
    endif

    do ns = 1,grid_size(1)
        dc%grad_A(ns, species) = (dc%A(ns+1, species)-dc%A(ns, species))/(g%s_vertices(ns+1)-g%s_vertices(ns))
    enddo

end subroutine calc_convection_coefficient_from_electric_field

function vertex_electric_field(i) result(electric_field)

    !-dPhi/ds [statV] at flux-surface vertex i: mean of the (1D) gradients of the two neighbouring shells,
    !the one shell at the axis and at the edge

    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_settings_mod, only: grid_size

    integer, intent(in) :: i
    real(dp) :: electric_field

    if (i.eq.1) then
        electric_field = -tetra_physics(1)%gPhi(1)
    elseif (i.eq.grid_size(1)+1) then
        electric_field = -tetra_physics((grid_size(1)-1)*6*grid_size(2)+1)%gPhi(1)
    else
        electric_field = -0.5_dp*(tetra_physics((i-2)*6*grid_size(2)+1)%gPhi(1) + &
                                  tetra_physics((i-1)*6*grid_size(2)+1)%gPhi(1))
    endif

end function vertex_electric_field

subroutine read_thermodiffusion_table()

    !Thermodiffusion coefficient c = D12/D11 - 3/2 of the electrons as a table in s and in the field -dPhi/ds [statV]
    !(GORILLA units, so no minor radius is needed). File format, '#' lines are comments:
    !   n_s n_E
    !   n_s*n_E lines "s  -dPhi/ds  c", s ascending in blocks of n_E lines, -dPhi/ds ascending within each block
    !Outside the table c is continued with the value at its edge.

    use gorilla_applets_types_mod, only: dc, in

    integer :: id, n_s, n_E, i, j, ios
    character(len=512) :: line
    real(dp) :: s_read

    open(newunit=id, file=trim(in%thermodiffusion_table_file), action='read', status='old')
    do
        read(id,'(a)') line
        line = adjustl(line)
        if (len_trim(line).eq.0) cycle
        if (line(1:1).ne.'#') exit
    enddo
    read(line,*) n_s, n_E
    if (allocated(dc%tab_s)) deallocate(dc%tab_s, dc%tab_E, dc%tab_c)
    allocate(dc%tab_s(n_s), dc%tab_E(n_E, n_s), dc%tab_c(n_E, n_s))
    do i = 1,n_s
        do j = 1,n_E
            do
                read(id,'(a)',iostat=ios) line
                if (ios.ne.0) then
                    print*, 'Error: thermodiffusion table ', trim(in%thermodiffusion_table_file), ' is too short'
                    stop
                endif
                line = adjustl(line)
                if ((len_trim(line).gt.0).and.(line(1:1).ne.'#')) exit
            enddo
            read(line,*) s_read, dc%tab_E(j,i), dc%tab_c(j,i)
            if (j.eq.1) then
                dc%tab_s(i) = s_read
            elseif (s_read.ne.dc%tab_s(i)) then
                print*, 'Error: thermodiffusion table: s changes inside a block (row ', i, ', column ', j, ')'
                stop
            endif
        enddo
    enddo
    close(id)
    dc%boole_thermo_table = .true.
    print*, 'thermodiffusion table read: ', n_s, ' values of s in [', dc%tab_s(1), ',', dc%tab_s(n_s), '], ', n_E, ' fields'

end subroutine read_thermodiffusion_table

function thermodiffusion_table_value(s0, field) result(c)

    !Bilinear interpolation of the thermodiffusion table (linear in the field within the two neighbouring rows,
    !then linear in s), constant beyond its edges

    use gorilla_applets_types_mod, only: dc

    real(dp), intent(in) :: s0, field
    real(dp) :: c, c_lo, c_hi, w
    integer :: i, n_s

    n_s = size(dc%tab_s)
    if (s0.le.dc%tab_s(1)) then
        c = interp_row(1)
    elseif (s0.ge.dc%tab_s(n_s)) then
        c = interp_row(n_s)
    else
        i = 1
        do while (dc%tab_s(i+1).lt.s0)
            i = i + 1
        enddo
        w = (s0 - dc%tab_s(i))/(dc%tab_s(i+1) - dc%tab_s(i))
        c_lo = interp_row(i)
        c_hi = interp_row(i+1)
        c = (1.0_dp - w)*c_lo + w*c_hi
    endif

contains

    function interp_row(ir) result(cr)
        integer, intent(in) :: ir
        real(dp) :: cr, wr
        integer :: j, n_E
        n_E = size(dc%tab_E, 1)
        if (field.le.dc%tab_E(1,ir)) then
            cr = dc%tab_c(1,ir)
        elseif (field.ge.dc%tab_E(n_E,ir)) then
            cr = dc%tab_c(n_E,ir)
        else
            j = 1
            do while (dc%tab_E(j+1,ir).lt.field)
                j = j + 1
            enddo
            wr = (field - dc%tab_E(j,ir))/(dc%tab_E(j+1,ir) - dc%tab_E(j,ir))
            cr = (1.0_dp - wr)*dc%tab_c(j,ir) + wr*dc%tab_c(j+1,ir)
        endif
    end function interp_row

end function thermodiffusion_table_value

subroutine sort_array(array)
  real(dp), dimension(:) :: array
  integer, parameter :: n = 6
  real :: a(n) = [3.5, 1.2, -4.0, 7.8, 0.0, 2.1]
  integer :: i, j
  real(dp) :: temp

  ! simple bubble sort (ascending)
  do i = 1, size(array)
     do j = i+1, size(array)
        if (array(j) < array(i)) then
           temp = array(i)
           array(i) = array(j)
           array(j) = temp
        end if
     end do
  end do

end subroutine sort_array

recursive subroutine quicksort(arr, left, right)
real(dp), intent(inout) :: arr(:)
integer, intent(in) :: left, right
integer :: i, j
real(dp) :: pivot, temp

i = left
j = right
pivot = arr((left+right)/2)

do
    do while (arr(i) < pivot)
        i = i + 1
    end do
    do while (arr(j) > pivot)
        j = j - 1
    end do
    if (i <= j) then
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
        i = i + 1
        j = j - 1
    end if
    if (i > j) exit
end do

if (left < j) call quicksort(arr, left, j)
if (i < right) call quicksort(arr, i, right)
end subroutine quicksort

end module utils_scef_diffusion_mod
