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
subroutine parallelised_particle_pushing_anomalous_transport(species)
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

                    call carry_out_anomalous_transport_displacement(x, ind_tetr, iface, t%step)
                endif

                ! Perform guiding-center orbit integration
                call orbit_timestep_anomalous_transport(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
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

end subroutine parallelised_particle_pushing_anomalous_transport

! ====================================================================
subroutine orbit_timestep_anomalous_transport(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
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

end subroutine orbit_timestep_anomalous_transport

! ====================================================================
subroutine carry_out_anomalous_transport_displacement(x, ind_tetr, iface, dt)
!
! Computes the displacement vector for anomalous transport and applies it.
!
! This subroutine calculates the perpendicular diffusion tensor and then
! calls displace_by_straight_line to perform the actual displacement
! through the tetrahedral mesh.
!
! The displacement follows:
!   Delta x^i = sqrt(2 * dt) * alpha^{ij} * xi_j + V_c^i * dt
!
! where x^1 = R, x^2 = phi, x^3 = Z (cylindrical coordinates),
! alpha is the Cholesky decomposition of D_perp, xi_j are random numbers
! with zero mean and unit variance, and V_c is the correction velocity.
!
! Input/Output:
!   x(3)      - Particle position in (R, phi, Z) coordinates
!   ind_tetr  - Current tetrahedron index (updated after displacement)
!   iface     - Current face index (updated after displacement)
!   dt        - Time step
!
    use tetra_physics_mod, only: tetra_physics
    use collis_ions, only: getran

    real(dp), dimension(3), intent(inout) :: x
    integer, intent(inout)                :: ind_tetr, iface
    real(dp), intent(in)                  :: dt

    real(dp), dimension(3) :: displacement, xi, V_c
    real(dp), parameter :: diffusion_coefficient = 1.0d-4  ! TODO: make this an input parameter

    ! Cholesky decomposition of the perpendicular diffusion tensor
    real(dp), dimension(3,3) :: alpha_perp_mat
    real(dp), dimension(3) :: h_contra, x_local
    real(dp) :: R_local, sqrt_2dt

    ! Compute position relative to first vertex
    x_local = x - tetra_physics(ind_tetr)%x1

    ! Interpolate R (major radius) at current position
    R_local = tetra_physics(ind_tetr)%R1 + dot_product(tetra_physics(ind_tetr)%gR, x_local)

    ! Interpolate h (unit vector along B) at current position using linear approximation
    ! h_i(x) = h_i_1 + grad(h_i) . (x - x1)
    h_contra(1) =  tetra_physics(ind_tetr)%h1_1 + dot_product(tetra_physics(ind_tetr)%gh1, x_local)
    h_contra(2) = (tetra_physics(ind_tetr)%h2_1 + dot_product(tetra_physics(ind_tetr)%gh2, x_local)) / (R_local**2)
    h_contra(3) =  tetra_physics(ind_tetr)%h3_1 + dot_product(tetra_physics(ind_tetr)%gh3, x_local)

    ! Compute Cholesky decomposition of the diffusion tensor
    call compute_diffusion_cholesky(h_contra, R_local, diffusion_coefficient, alpha_perp_mat)

    ! Compute correction velocity
    call compute_correction_velocity(ind_tetr, h_contra, R_local, diffusion_coefficient, V_c)

    ! Generate random numbers with zero mean and unit variance
    call getran(0, xi(1))
    call getran(0, xi(2))
    call getran(0, xi(3))

    ! Compute displacement: Delta x^i = sqrt(2 * dt) * alpha^{ij} * xi_j + V_c^i * dt
    ! Note: alpha_perp_mat is lower triangular, so alpha^{ij} * xi_j expands as:
    ! displacement(1) = alpha(1,1)*xi(1)
    ! displacement(2) = alpha(2,1)*xi(1) + alpha(2,2)*xi(2)
    ! displacement(3) = alpha(3,1)*xi(1) + alpha(3,2)*xi(2) + alpha(3,3)*xi(3)
    sqrt_2dt = sqrt(2.0_dp * dt)
    displacement(1) = sqrt_2dt * alpha_perp_mat(1,1) * xi(1) + V_c(1) * dt
    displacement(2) = sqrt_2dt * (alpha_perp_mat(2,1) * xi(1) + alpha_perp_mat(2,2) * xi(2)) + V_c(2) * dt
    displacement(3) = sqrt_2dt * (alpha_perp_mat(3,1) * xi(1) + alpha_perp_mat(3,2) * xi(2) &
                                + alpha_perp_mat(3,3) * xi(3)) + V_c(3) * dt

    call displace_by_straight_line(x, ind_tetr, iface, displacement)

end subroutine carry_out_anomalous_transport_displacement

! ====================================================================
subroutine compute_diffusion_cholesky(h_contra, R_local, D_perp, alpha_perp_mat)
!
! Computes the Cholesky decomposition of the perpendicular diffusion tensor.
!
! The perpendicular diffusion tensor is:
!   D_perp^{ij} = D_perp * (g^{ij} - h^i * h^j)
!
! where g^{ij} is the contravariant metric tensor (cylindrical coordinates:
! g^{11} = g^{33} = 1, g^{22} = 1/R^2, off-diagonal = 0) and h^i are the
! contravariant components of the unit vector in the magnetic field direction.
!
! The Cholesky decomposition produces a lower triangular matrix alpha such that
! D_perp = alpha * alpha^T.
!
! Input:
!   h_contra(3) - Contravariant components of the magnetic field unit vector
!   R_local     - Local major radius
!   D_perp      - Perpendicular diffusion coefficient (scalar)
!
! Output:
!   alpha_perp_mat(3,3) - Lower triangular Cholesky factor
!
    real(dp), dimension(3), intent(in)    :: h_contra
    real(dp), intent(in)                  :: R_local, D_perp
    real(dp), dimension(3,3), intent(out) :: alpha_perp_mat

    real(dp), dimension(3,3) :: D_perp_mat, g_mat

    ! Contravariant metric tensor for cylindrical coordinates (R, phi, Z)
    ! g^{11} = 1, g^{22} = 1/R^2, g^{33} = 1, off-diagonal = 0
    g_mat = 0.0_dp
    g_mat(1,1) = 1.0_dp
    g_mat(2,2) = 1.0_dp / (R_local**2)
    g_mat(3,3) = 1.0_dp

    ! Perpendicular diffusion tensor: D_perp^{ij} = D_perp * (g^{ij} - h^i * h^j)
    D_perp_mat(1,1) = D_perp * (g_mat(1,1) - h_contra(1) * h_contra(1))
    D_perp_mat(1,2) = D_perp * (g_mat(1,2) - h_contra(1) * h_contra(2))
    D_perp_mat(1,3) = D_perp * (g_mat(1,3) - h_contra(1) * h_contra(3))
    D_perp_mat(2,1) = D_perp * (g_mat(2,1) - h_contra(2) * h_contra(1))
    D_perp_mat(2,2) = D_perp * (g_mat(2,2) - h_contra(2) * h_contra(2))
    D_perp_mat(2,3) = D_perp * (g_mat(2,3) - h_contra(2) * h_contra(3))
    D_perp_mat(3,1) = D_perp * (g_mat(3,1) - h_contra(3) * h_contra(1))
    D_perp_mat(3,2) = D_perp * (g_mat(3,2) - h_contra(3) * h_contra(2))
    D_perp_mat(3,3) = D_perp * (g_mat(3,3) - h_contra(3) * h_contra(3))

    ! Cholesky decomposition: alpha_perp such that D_perp = alpha_perp * alpha_perp^T
    ! Row 1
    alpha_perp_mat(1,1) = sqrt(D_perp_mat(1,1))
    alpha_perp_mat(1,2) = 0.0_dp
    alpha_perp_mat(1,3) = 0.0_dp
    ! Row 2
    alpha_perp_mat(2,1) = D_perp_mat(1,2) / alpha_perp_mat(1,1)
    alpha_perp_mat(2,2) = sqrt(D_perp_mat(2,2) - alpha_perp_mat(2,1)**2)
    alpha_perp_mat(2,3) = 0.0_dp
    ! Row 3
    alpha_perp_mat(3,1) = D_perp_mat(1,3) / alpha_perp_mat(1,1)
    alpha_perp_mat(3,2) = (D_perp_mat(2,3) - alpha_perp_mat(2,1) * alpha_perp_mat(3,1)) / alpha_perp_mat(2,2)
    alpha_perp_mat(3,3) = sqrt(D_perp_mat(3,3) - alpha_perp_mat(3,1)**2 - alpha_perp_mat(3,2)**2)

end subroutine compute_diffusion_cholesky

! ====================================================================
subroutine compute_correction_velocity(ind_tetr, h_contra, R_local, D_perp, V_c)
!
! Computes the correction velocity for the stochastic differential equation.
!
! The correction velocity is:
!   V_c^i = V^i + (1/sqrt(g)) * d/dx^j (sqrt(g) * D_perp^{ij})
!         = V^i + dD_perp^{ij}/dx^j + D_perp^{i1}/R
!
! where V^i is the drift velocity (set to zero for now), and the second term
! comes from the coordinate system (cylindrical: sqrt(g) = R).
!
! The derivative of the diffusion tensor is:
!   dD_perp^{ij}/dx^j = -D_perp * (h^i * dh^j/dx^j + h^j * dh^i/dx^j)
!
! Input:
!   ind_tetr    - Current tetrahedron index
!   h_contra(3) - Contravariant components of the magnetic field unit vector
!   R_local     - Local major radius
!   D_perp      - Perpendicular diffusion coefficient (scalar)
!
! Output:
!   V_c(3)      - Correction velocity
!
    use tetra_physics_mod, only: tetra_physics

    integer, intent(in)                   :: ind_tetr
    real(dp), dimension(3), intent(in)    :: h_contra
    real(dp), intent(in)                  :: R_local, D_perp
    real(dp), dimension(3), intent(out)   :: V_c

    real(dp), dimension(3) :: V, partial_D_perp_ij_partial_xj, D_perp_i1
    real(dp), dimension(3,3) :: grad_h_contra  ! grad_h_contra(i,j) = dh^i/dx^j
    real(dp) :: div_h_contra

    ! Drift velocity (set to zero for now)
    V = 0.0_dp

    ! Gradient of contravariant h components: grad_h_contra(i,j) = dh^i/dx^j
    ! For cylindrical coordinates, h^i = g^{ij} h_j, so:
    ! h^1 = h_1, h^2 = h_2/R^2, h^3 = h_3
    ! dh^1/dx^j = dh_1/dx^j = gh1(j)
    ! dh^2/dx^j = d(h_2/R^2)/dx^j = gh2(j)/R^2 - 2*h^2/R * dR/dx^j
    ! dh^3/dx^j = dh_3/dx^j = gh3(j)
    grad_h_contra(1,:) = tetra_physics(ind_tetr)%gh1
    grad_h_contra(2,:) = tetra_physics(ind_tetr)%gh2 / (R_local**2) &
                       - 2.0_dp * h_contra(2) / R_local * tetra_physics(ind_tetr)%gR
    grad_h_contra(3,:) = tetra_physics(ind_tetr)%gh3

    ! Divergence of h^j: div_h_contra = dh^j/dx^j (sum over j)
    div_h_contra = grad_h_contra(1,1) + grad_h_contra(2,2) + grad_h_contra(3,3)

    ! Derivative of perpendicular diffusion tensor: dD_perp^{ij}/dx^j = -D_perp * (h^i * dh^j/dx^j + h^j * dh^i/dx^j)
    ! Summed over j, this gives a vector (one component for each i)
    partial_D_perp_ij_partial_xj(1) = -D_perp * ( &
        h_contra(1) * div_h_contra + &
        h_contra(1) * grad_h_contra(1,1) + h_contra(2) * grad_h_contra(1,2) + h_contra(3) * grad_h_contra(1,3))
    partial_D_perp_ij_partial_xj(2) = -D_perp * ( &
        h_contra(2) * div_h_contra + &
        h_contra(1) * grad_h_contra(2,1) + h_contra(2) * grad_h_contra(2,2) + h_contra(3) * grad_h_contra(2,3))
    partial_D_perp_ij_partial_xj(3) = -D_perp * ( &
        h_contra(3) * div_h_contra + &
        h_contra(1) * grad_h_contra(3,1) + h_contra(2) * grad_h_contra(3,2) + h_contra(3) * grad_h_contra(3,3))

    ! Compute D_perp^{i1} = D_perp * (g^{i1} - h^i * h^1)
    ! For cylindrical coordinates: g^{11} = 1, g^{21} = g^{31} = 0
    ! So: D_perp^{11} = D_perp * (1 - h^1 * h^1)
    !     D_perp^{21} = D_perp * (0 - h^2 * h^1) = -D_perp * h^2 * h^1
    !     D_perp^{31} = D_perp * (0 - h^3 * h^1) = -D_perp * h^3 * h^1
    D_perp_i1(1) = D_perp * (1.0_dp - h_contra(1) * h_contra(1))
    D_perp_i1(2) = -D_perp * h_contra(2) * h_contra(1)
    D_perp_i1(3) = -D_perp * h_contra(3) * h_contra(1)

    ! Correction velocity: V_c^i = V^i + dD_perp^{ij}/dx^j + D_perp^{i1}/R
    V_c(1) = V(1) + partial_D_perp_ij_partial_xj(1) + D_perp_i1(1) / R_local
    V_c(2) = V(2) + partial_D_perp_ij_partial_xj(2) + D_perp_i1(2) / R_local
    V_c(3) = V(3) + partial_D_perp_ij_partial_xj(3) + D_perp_i1(3) / R_local

end subroutine compute_correction_velocity

! ====================================================================
subroutine displace_by_straight_line(x, ind_tetr, iface, displacement)
!
! Displaces a particle along a straight line through the tetrahedral mesh.
!
! This version traverses tetrahedra along the displacement path instead of
! using find_tetra after each displacement. This is much more efficient as
! it only needs to check faces of the current tetrahedron rather than
! searching through all tetrahedra.
!
! Input/Output:
!   x(3)        - Particle position in (s, theta, phi) coordinates
!   ind_tetr    - Current tetrahedron index (updated after displacement)
!   iface       - Current face index (updated after displacement)
!   displacement(3) - Displacement vector in (s, theta, phi) coordinates
!
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: tetra_grid
    use pusher_tetra_func_mod, only: pusher_handover2neighbour
    use tetra_grid_settings_mod, only: sfc_s_min
    use find_tetra_mod, only: find_tetra
    use constants, only: eps

    real(dp), dimension(3), intent(inout) :: x
    integer, intent(inout)                :: ind_tetr, iface
    real(dp), dimension(3), intent(in)    :: displacement

    real(dp) :: remaining_distance, distance_to_face, dist_min, total_distance
    real(dp), dimension(3) :: direction, z_rel
    integer :: hit_face, iper_phi, ind_tetr_old
    integer :: max_crossings, crossing_count
    real(dp) :: vpar_dummy, vperp_dummy

    ! Compute total displacement distance and direction
    total_distance = sqrt(displacement(1)**2 + displacement(2)**2 + displacement(3)**2)
    if (total_distance < 1.0d-15) return  ! No displacement

    direction = displacement / total_distance
    remaining_distance = total_distance

    ! Maximum number of tetrahedron crossings to prevent infinite loops
    max_crossings = 10000
    crossing_count = 0

    ! Use relative tolerance based on tetrahedron size, consistent with find_tetra_face_intersection
    dist_min = eps * abs(tetra_physics(ind_tetr)%dist_ref)

    ! Traverse tetrahedra until displacement is complete
    do while (remaining_distance > dist_min .and. crossing_count < max_crossings)
        crossing_count = crossing_count + 1

        ! Compute position relative to first vertex of current tetrahedron
        z_rel = x - tetra_physics(ind_tetr)%x1

        ! Find intersection with tetrahedron faces
        call find_tetra_face_intersection(z_rel, direction, ind_tetr, distance_to_face, hit_face)

        if (hit_face == 0) then
            ! No valid intersection found - this shouldn't happen if particle is inside
            ! Move the full remaining distance and exit
            x = x + remaining_distance * direction
            exit
        endif

        if (distance_to_face >= remaining_distance) then
            ! Displacement ends within this tetrahedron
            x = x + remaining_distance * direction
            remaining_distance = 0.0_dp
        else
            ! Move to the face and continue to next tetrahedron
            x = x + distance_to_face * direction

            ! Get the neighboring tetrahedron
            ind_tetr_old = ind_tetr
            call pusher_handover2neighbour(ind_tetr_old, ind_tetr, hit_face, x, iper_phi)

            ! Update iface to the entry face of the new tetrahedron
            iface = hit_face

            ! Check if we've left the domain
            if (ind_tetr == -1) then
                ! Handle boundary: reflect at inner boundary, mark lost at outer
                if (x(1) < 1.01_dp * sfc_s_min) then
                    ! Reflect at inner boundary
                    x(1) = 2.0_dp * sfc_s_min - x(1)
                    direction(1) = -direction(1)  ! Reverse radial direction
                    ! Need to find tetrahedron again after reflection
                    vpar_dummy = 0.0_dp
                    vperp_dummy = 0.0_dp
                    call find_tetra(x, vpar_dummy, vperp_dummy, ind_tetr, iface)
                    if (ind_tetr == -1) exit  ! Lost particle
                else
                    ! Lost at outer boundary
                    exit
                endif
            endif

            remaining_distance = remaining_distance - distance_to_face
        endif
    enddo

    ! Final boundary checks
    if (x(1) < sfc_s_min) then
        x(1) = 2.0_dp * sfc_s_min - x(1)
    endif
    if (x(1) > 1.0_dp) then
        x(1) = 1.0_dp
        ind_tetr = -1  ! Mark as lost
    endif

end subroutine displace_by_straight_line

! ====================================================================
subroutine find_tetra_face_intersection(start_point, direction, ind_tetr, &
                                        distance_to_face, hit_face)
!
! Finds the intersection of a ray with a tetrahedron's faces.
!
! Uses the pre-computed face normals from tetra_physics to find which face
! the ray exits through. This is analogous to how the particle pusher works.
!
! The tetrahedron face normals point INWARD. A particle inside the tetrahedron
! will cross face i when the normal distance to face i becomes negative.
! For a ray from start_point in direction 'direction', we find when:
!   normal_distance(t) = normal_distance_0 + t * (anorm . direction) = 0
!
! Special handling for particles on a face:
!   When a particle is exactly on a face (normal_dist_0 ≈ 0), we need to
!   determine if it should exit through that face. We accept t_intersect >= 0
!   but when multiple faces have t ≈ 0, we pick the one we're moving towards
!   most strongly (largest negative normal_dot_dir).
!
! Face indexing in GORILLA:
!   Face 1: opposite to vertex 1 (contains vertices 2,3,4)
!   Face 2: opposite to vertex 2 (contains vertices 1,3,4)
!   Face 3: opposite to vertex 3 (contains vertices 1,2,4)
!   Face 4: opposite to vertex 4 (contains vertices 1,2,3)
!
! Input:
!   start_point(3)    - Starting position inside the tetrahedron (relative to x1)
!   direction(3)      - Direction vector of the ray (NOT necessarily unit)
!   ind_tetr          - Index of the current tetrahedron
!
! Output:
!   distance_to_face  - Distance along direction to the exit face
!                       (in units of |direction|, so actual distance = distance_to_face * |direction|)
!   hit_face          - Index of the exit face (1-4), or 0 if no valid intersection
!
    use tetra_physics_mod, only: tetra_physics
    use constants, only: eps

    real(dp), intent(in)  :: start_point(3)
    real(dp), intent(in)  :: direction(3)
    integer,  intent(in)  :: ind_tetr
    real(dp), intent(out) :: distance_to_face
    integer,  intent(out) :: hit_face

    real(dp) :: normal_dot_dir, normal_dist_0, t_intersect, dist_min
    real(dp) :: best_normal_dot_dir
    integer  :: iface

    hit_face = 0
    distance_to_face = huge(distance_to_face)
    best_normal_dot_dir = 0.0_dp

    ! Use relative tolerance based on tetrahedron size, consistent with isinside()
    dist_min = eps * abs(tetra_physics(ind_tetr)%dist_ref)

    ! Loop over all 4 faces
    do iface = 1, 4
        ! Compute: anorm . direction
        ! If negative, the ray is moving towards this face (since normals point inward)
        normal_dot_dir = dot_product(tetra_physics(ind_tetr)%anorm(:,iface), direction)

        ! Only consider faces we're moving towards (normal_dot_dir < 0 means moving towards the face)
        if (normal_dot_dir >= -dist_min) cycle

        ! Compute current normal distance to this face
        ! For face 1: includes dist_ref offset
        ! For faces 2-4: just anorm . x_local where x_local = position relative to x1
        if (iface == 1) then
            normal_dist_0 = dot_product(tetra_physics(ind_tetr)%anorm(:,1), start_point) &
                          + tetra_physics(ind_tetr)%dist_ref
        else
            normal_dist_0 = dot_product(tetra_physics(ind_tetr)%anorm(:,iface), start_point)
        endif

        ! Solve for t: normal_dist_0 + t * normal_dot_dir = 0
        ! => t = -normal_dist_0 / normal_dot_dir
        t_intersect = -normal_dist_0 / normal_dot_dir

        ! Accept intersections in forward direction (t >= -dist_min to handle numerical noise)
        if (t_intersect < -dist_min) cycle

        ! Determine if this is the best candidate
        if (t_intersect < distance_to_face - dist_min) then
            ! Clearly closer than current best
            distance_to_face = t_intersect
            hit_face = iface
            best_normal_dot_dir = normal_dot_dir
        else if (abs(t_intersect - distance_to_face) <= dist_min) then
            ! Similar distance (particle near an edge/corner) - pick face we're moving towards most strongly
            if (normal_dot_dir < best_normal_dot_dir) then
                distance_to_face = t_intersect
                hit_face = iface
                best_normal_dot_dir = normal_dot_dir
            endif
        endif
    enddo

end subroutine find_tetra_face_intersection

end module utils_anomalous_transport_mod
