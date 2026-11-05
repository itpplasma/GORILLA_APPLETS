module utils_scef_diffusion_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

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
        'A_fit_first_half A_fit_second_half rms_residual_of_<ds>(t)_fit B_fit B_end sigma_B_end n_passes'
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
        write(quality_unit,'(i4,25es14.5,i3)') ns, s%s0, start%t(species), n_coll, dble(s%n_particles), &
            dble(n_lost)/s%n_particles, mean_x - s%s0, std_x, sqrt(2.0_dp*max(B,0.0_dp)*start%t(species)), std_x/cell_size, &
            skewness, kurtosis, q, A, A_end, sigma_A_end, A_half1, A_half2, rms_residual, B, B_end, sigma_B_end, pass
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
        write(file_id,*) dc%A(i, species), dc%B(i, species)
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

    real(dp) :: delta_x, delta_t, xi, A, B, cell_size, s_c, B_c, s_lo, s_hi, c_thermo
    integer :: i, ns, k, num_steps_min, count_lost_particles, num_particles, particle_multiplication
    integer :: density_unit, n, ind_tetr, iface
    real(dp), dimension(3) :: x_particle, z_save
    real(dp) :: vpar_init, vperp_init
    real(dp), dimension(:), allocatable :: rw_density, rw_prism_densities
    real(dp), dimension(:), allocatable :: position, weight, exit_time
    real(dp), dimension(:), allocatable :: marker_shell   !density contribution per flux shell of one marker (mean over its walkers)
    type(time_t) :: t
    logical :: boole_lost, boole_use_fit_function = .true.

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


    if (.not.dc%boole_cached_data_loaded(species)) then
        call get_diffusion_coefficient_data(species, boole_use_fit_function)
        dc%boole_cached_data_loaded(species) = .true.
    endif
    call calc_convection_coefficient_from_electric_field(species, boole_use_fit_function)

    t%step = in%time_step
    num_steps_min = 1000

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
        if (mod(i-1,particle_multiplication).eq.0) marker_shell = 0.0_dp
        do while ((t%confined.lt.t%step).and.(.not.boole_lost))

            !binsrc finds k such that g%s_vertices(k-1) < position(i) < g%s_vertices(k)
            call binsrc(g%s_vertices,1,grid_size(1)+1,position(i),k)
            k = k-1
            cell_size = g%s_vertices(k+1) - g%s_vertices(k)
            A = dc%A(k, species) + dc%grad_A(k, species)*(position(i)-g%s_vertices(k))
            B = dc%B(k, species) + dc%grad_B(k, species)*(position(i)-g%s_vertices(k))

            if (boole_use_fit_function) then
                !B from the fit of ln B; A holds q*E/T here, so the E_r drift uses the same B as the
                !diffusion. Below s_c, B -> 0 like s (B = (ds/dr)^2 D_r with finite D_r). Since the
                !volume element is constant in s, the walk needs the Ito drift dB/ds.
                s_c = g%s_vertices(4)
                if (position(i).lt.s_c) then
                    B_c = exp(sum(dc%polynomial_coefficients_for_B(:, species)*(/1.0_dp,s_c,s_c**2,0.0_dp/)))
                    B = B_c*position(i)/s_c
                    A = B*A + B_c/s_c
                else
                    B = exp(sum(dc%polynomial_coefficients_for_B(:, species) &
                                *(/1.0_dp,position(i),position(i)**2,0.0_dp/)))
                    A = B*A + B*sum(dc%polynomial_coefficients_for_A(:, species) &
                                    *(/1.0_dp,position(i),0.0_dp/)) !dB/ds = B * d(ln B)/ds
                endif
            endif

            !Thermodiffusion: the walkers carry no energy, so the flux driven by the temperature gradient,
            !-c*B*n*d(ln T)/ds with c = D12/D11 - 3/2, enters as an explicit drift (central difference in s)
            !With thermodiffusion_table_file, c comes from the table at the current field (linear between vertices)
            c_thermo = in%thermodiffusion_coefficient
            if (dc%boole_thermo_table.and.(species.eq.2)) c_thermo = dc%c_thermo(k, species) &
                + (dc%c_thermo(k+1, species)-dc%c_thermo(k, species))*(position(i)-g%s_vertices(k))/cell_size
            if (c_thermo.ne.0.0_dp) then
                s_lo = max(position(i) - 1.0d-4, 0.0_dp)
                s_hi = min(position(i) + 1.0d-4, 1.0_dp)
                A = A - c_thermo*B &
                        *log(target_temperature(rad_from_s(s_hi), species)/target_temperature(rad_from_s(s_lo), species)) &
                        /(s_hi - s_lo)
            endif


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
            if (position(i).gt.1.0_dp)    then
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

end subroutine calc_density_via_random_walk

subroutine get_diffusion_coefficient_data(species, boole_use_fit_function)

    use gorilla_applets_types_mod, only: dc, g, in
    use tetra_grid_settings_mod, only: grid_size
    use utils_polyfit_mod, only: quadratic_fit

    integer, intent(in) :: species
    logical, intent(in) :: boole_use_fit_function

    character(len=100) :: filename, species_str
    integer :: id, ns, lower_bound, upper_bound
    real(dp), dimension(:), allocatable :: x,y
    logical :: success
    real(dp) :: a0, a1, a2


    if (.not.allocated(dc%A))                allocate(dc%A(grid_size(1)+1, 2))
    if (.not.allocated(dc%A_from_first_run)) allocate(dc%A_from_first_run(grid_size(1)+1, 2))
    if (.not.allocated(dc%B))                allocate(dc%B(grid_size(1)+1, 2))
    if (.not.allocated(dc%grad_A))           allocate(dc%grad_A(grid_size(1), 2))
    if (.not.allocated(dc%grad_B))           allocate(dc%grad_B(grid_size(1), 2))

    write(species_str, '(I0)') species
    filename = 'A_and_B_species' // trim(species_str) // '.dat'
    open(newunit=id,file = filename,action='read')
    do ns = 1,grid_size(1)+1
        read(id,*) dc%A(ns, species), dc%B(ns, species)
    enddo
    close(id)

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
        if (i.eq.1) then
            electric_field = -tetra_physics(1)%gPhi(1)
        elseif (i.eq.grid_size(1)+1) then
            electric_field = -tetra_physics((grid_size(1)-1)*6*grid_size(2)+1)%gPhi(1)
        else
            electric_field = -0.5_dp*(tetra_physics((i-2)*6*grid_size(2)+1)%gPhi(1) + &
                                      tetra_physics((i-1)*6*grid_size(2)+1)%gPhi(1))
        endif

        !q*E/T; with the fit function the random walk multiplies it by its fitted B, otherwise by the
        !tabulated B here
        dc%A(i, species) = start%particle_charge(species)*electric_field &
                           /(target_temperature(rad_from_s(g%s_vertices(i)), species)*ev2erg)
        if (.not.boole_use_fit_function) dc%A(i, species) = dc%B(i, species)*dc%A(i, species) &
                                                            + dc%A_from_first_run(i, species)
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
