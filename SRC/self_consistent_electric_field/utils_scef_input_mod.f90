module utils_scef_input_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine read_self_consistent_electric_field_inp_into_type

    use gorilla_applets_types_mod, only: in

    real(dp) :: time_step,energy_eV,n_particles, density, source_relaxation_factor, quasineutrality_tol, diffusion_scan_s0, diffusion_tau, &
                diffusion_tau_collisions, edge_buffer_width, noise_tolerance_factor, potential_relaxation_factor
    character(len=256) :: restart_phi_file, restart_source_file
    logical :: boole_adaptive_relaxation, boole_measured_drift, boole_continuous_background
    integer :: source_buffer_mode, collision_background_mode
    real(dp) :: source_buffer_bound_factor
    integer :: recalibration_mode
    real(dp) :: thermodiffusion_coefficient, ion_tracing_time_limit, diffusion_sigma_A_target, diffusion_max_spread_cells
    character(len=256) :: thermodiffusion_table_file
    logical ::boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_monoenergetic, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
               boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
               boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, boole_static_ne, reduced_output
    integer :: i_integrator_type, seed_option, n_electric_potential_updates, update_dimension, n_species, n_source_updates
    integer :: i_collision_mode, n_particles_diffusion
    logical :: boole_honest_tracing(2), boole_recompute_D(2)

    integer :: s_inp_unit

    !Namelist for self consistent electric field input
    NAMELIST /self_consistent_ef_nml/ time_step,energy_eV,n_particles,boole_squared_moments,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_monoenergetic, boole_linear_density_simulation, &
    & boole_antithetic_variate,boole_linear_temperature_simulation,i_integrator_type,seed_option, boole_write_vertex_indices, &
    & boole_write_vertex_coordinates, boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
    & boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_electric_potential_updates, update_dimension, &
    & n_species, boole_static_ne, boole_honest_tracing, boole_recompute_D, n_source_updates, reduced_output, source_relaxation_factor, &
    & quasineutrality_tol, i_collision_mode, n_particles_diffusion, diffusion_scan_s0, diffusion_tau, &
    & diffusion_tau_collisions, edge_buffer_width, noise_tolerance_factor, potential_relaxation_factor, &
    & boole_adaptive_relaxation, restart_phi_file, restart_source_file, thermodiffusion_coefficient, &
    & ion_tracing_time_limit, diffusion_sigma_A_target, diffusion_max_spread_cells, thermodiffusion_table_file, &
    & boole_measured_drift, recalibration_mode, boole_continuous_background, source_buffer_mode, &
    & source_buffer_bound_factor, collision_background_mode

    quasineutrality_tol = 0.0_dp   !optional in the .inp: no early exit unless set
    i_collision_mode = 1           !optional: full collision operator
    n_particles_diffusion = 40000  !optional
    diffusion_scan_s0 = -1.0_dp    !optional: calibrate all flux surfaces, then run SCEF
    diffusion_tau = -1.0_dp        !optional: hardcoded calibration time per species
    diffusion_tau_collisions = -1.0_dp   !optional: calibration time from diffusion_tau
    edge_buffer_width = 0.0_dp     !optional: no outer buffer
    noise_tolerance_factor = 0.0_dp   !optional: no noise-based stopping
    potential_relaxation_factor = 1.0_dp/3.0_dp   !optional
    boole_adaptive_relaxation = .false.   !optional: one global relaxation factor
    thermodiffusion_coefficient = 0.0_dp  !optional: no thermodiffusive drift in the random walk
    ion_tracing_time_limit = -1.0_dp      !optional: ions traced up to time_step
    restart_phi_file = ''        !optional: start from the initial potential
    restart_source_file = ''     !optional: start from the target density as source
    diffusion_sigma_A_target = -1.0_dp    !optional: calibration time from the collision time only
    diffusion_max_spread_cells = 1.0_dp   !optional
    thermodiffusion_table_file = ''       !optional: constant thermodiffusion_coefficient
    boole_measured_drift = .false.        !optional: walk with the fitted ln B and the field drift only
    recalibration_mode = 0                !optional: 0 calibrate A,B once; 1 also before every later source step; 2 before every later potential update
    boole_continuous_background = .false.   !optional, legacy: .true. = collision_background_mode 1
    collision_background_mode = 0           !optional: 0 per prism (staircase), 1 exact at the particle, 2 linear in the tetrahedron
    source_buffer_mode = 0                !optional: 0 source in the edge buffer not updated; 1 updated, |S| capped
    source_buffer_bound_factor = 2.0_dp   !optional: cap for source_buffer_mode 1, times max |S| outside the buffer
    open(newunit = s_inp_unit, file='self_consistent_ef.inp', status='unknown')
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
    in%boole_honest_tracing = boole_honest_tracing
    in%boole_recompute_D    = boole_recompute_D
    in%n_source_updates         = n_source_updates
    in%reduced_output           = reduced_output
    in%source_relaxation_factor = source_relaxation_factor
    in%quasineutrality_tol      = quasineutrality_tol
    in%i_collision_mode         = i_collision_mode
    in%n_particles_diffusion    = n_particles_diffusion
    in%diffusion_scan_s0        = diffusion_scan_s0
    in%diffusion_tau            = diffusion_tau
    in%diffusion_tau_collisions = diffusion_tau_collisions
    in%edge_buffer_width        = edge_buffer_width
    in%noise_tolerance_factor   = noise_tolerance_factor
    in%potential_relaxation_factor = potential_relaxation_factor
    in%boole_adaptive_relaxation = boole_adaptive_relaxation
    in%thermodiffusion_coefficient = thermodiffusion_coefficient
    in%ion_tracing_time_limit = ion_tracing_time_limit
    in%restart_phi_file = restart_phi_file
    in%restart_source_file = restart_source_file
    in%diffusion_sigma_A_target = diffusion_sigma_A_target
    in%diffusion_max_spread_cells = diffusion_max_spread_cells
    in%thermodiffusion_table_file = thermodiffusion_table_file
    in%boole_measured_drift = boole_measured_drift
    in%recalibration_mode = recalibration_mode
    in%collision_background_mode = collision_background_mode
    if (boole_continuous_background .and. (collision_background_mode.eq.0)) in%collision_background_mode = 1
    in%source_buffer_mode = source_buffer_mode
    in%source_buffer_bound_factor = source_buffer_bound_factor

    print *,'GORILLA_APPLETS: Loaded input data from self_consistent_ef.inp'

end subroutine read_self_consistent_electric_field_inp_into_type

subroutine print_errors_for_bad_inputs

    use gorilla_applets_types_mod, only: in
    use tetra_physics_mod, only: coord_system

    if (in%i_integrator_type.eq.2) then
        print*, 'Error: i_integrator_type set to 2, this module only works with i_integrator_type set to 1'
        print*, 'Program terminated'
        stop
    endif

    if (in%boole_refined_sqrt_g.eqv..true.) then
        print*, 'Error: boole_refined_sqrt_g set to .true., but this only works for cylindrical coordinates. This module &
                 works with flux coordinates.'
        print*, 'Program terminated'
        stop
    endif

    if (in%boole_write_refined_prism_volumes.eqv..true.) then
        print*, 'Error: boole_write_refined_prism_volumes set to .true., but this only works for cylindrical coordinates. &
                 This module works with flux coordinates.'
        print*, 'Program terminated'
        stop
    endif

    if (coord_system.eq.1) then
        print*, 'Error: coord_system set to 1, but this module works with flux coordinates.'
        print*, 'Program terminated'
        stop
    endif

    if (in%n_source_updates.lt.0) then
        print*, 'Error: n_source_updates must be >= 0 (0 = fixed source, no source tuning).'
        print*, 'Program terminated'
        stop
    endif

    if ((in%update_dimension.ne.1) .and. (in%update_dimension.ne.3)) then
        print*, 'Error: update_dimension must be 1 (flux-surface updates) or 3 (3D updates).'
        print*, 'Program terminated'
        stop
    endif

    if ((in%update_dimension.eq.3) .and. (.not.all(in%boole_honest_tracing))) then
        print*, 'Error: update_dimension = 3 requires boole_honest_tracing = .true. for both species; &
                 &the random walk only resolves the radial direction.'
        print*, 'Program terminated'
        stop
    endif

    if ((in%update_dimension.eq.3) .and. (in%n_source_updates.gt.0)) then
        print*, 'Error: update_dimension = 3 is not implemented together with source tuning; set n_source_updates = 0.'
        print*, 'Program terminated'
        stop
    endif

    if ((in%update_dimension.eq.3) .and. in%boole_static_ne) then
        print*, 'Error: update_dimension = 3 is not implemented for boole_static_ne = .true.'
        print*, 'Program terminated'
        stop
    endif

    if (in%quasineutrality_tol.lt.0.0_dp) then
        print*, 'Error: quasineutrality_tol must be >= 0 (0 = no early exit).'
        print*, 'Program terminated'
        stop
    endif

    if (all(in%i_collision_mode.ne.[1,2,3,4])) then
        print*, 'Error: i_collision_mode must be 1 (full), 2 (energy scattering and drag), 3 (drag) or 4 (pitch-angle only).'
        print*, 'Program terminated'
        stop
    endif

    if ((in%diffusion_scan_s0.gt.0.0_dp) .and. &
        (.not.any((.not.in%boole_honest_tracing).and.in%boole_recompute_D))) then
        print*, 'Error: diffusion_scan_s0 > 0 needs a species with boole_honest_tracing = .false. and boole_recompute_D = .true.'
        print*, 'Program terminated'
        stop
    endif

    if ((in%potential_relaxation_factor.le.0.0_dp) .or. (in%potential_relaxation_factor.gt.1.0_dp)) then
        print*, 'Error: potential_relaxation_factor must be in (0, 1].'
        print*, 'Program terminated'
        stop
    endif

    if ((in%edge_buffer_width.lt.0.0_dp) .or. (in%edge_buffer_width.ge.1.0_dp)) then
        print*, 'Error: edge_buffer_width must be in [0, 1) (width in r/a, 0 = off).'
        print*, 'Program terminated'
        stop
    endif

    if ((in%noise_tolerance_factor.lt.0.0_dp) .or. &
        ((in%noise_tolerance_factor.gt.0.0_dp) .and. (in%update_dimension.ne.1))) then
        print*, 'Error: noise_tolerance_factor must be >= 0, and > 0 is only implemented for update_dimension = 1.'
        print*, 'Program terminated'
        stop
    endif

    if ((in%edge_buffer_width.gt.0.0_dp) .and. (in%update_dimension.ne.1)) then
        print*, 'Error: edge_buffer_width > 0 is only implemented for update_dimension = 1.'
        print*, 'Program terminated'
        stop
    endif

    if ((in%boole_adaptive_relaxation .or. (len_trim(in%restart_phi_file).gt.0)) .and. &
        ((in%update_dimension.ne.1) .or. in%boole_static_ne)) then
        print*, 'Error: boole_adaptive_relaxation and restart_phi_file need update_dimension = 1 and boole_static_ne = .false.'
        print*, 'Program terminated'
        stop
    endif

    if ((in%collision_background_mode.lt.0) .or. (in%collision_background_mode.gt.2)) then
        print*, 'Error: collision_background_mode must be 0, 1 or 2'
        stop
    endif

    if ((in%collision_background_mode.gt.0) .and. ((coord_system.ne.2) .or. &
        in%boole_preserve_energy_and_momentum_during_collisions)) then
        print*, 'Error: collision_background_mode 1 and 2 (and boole_continuous_background) need coord_system = 2 &
                &(s = x(1)) and boole_preserve_energy_and_momentum_during_collisions = .false.'
        stop
    endif

end subroutine print_errors_for_bad_inputs

end module utils_scef_input_mod
