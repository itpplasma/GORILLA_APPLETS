module utils_rmp_response_currents_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64
    ! Canonical KIM r_eff-alignment flag lives in profile_data_mod; use-associate
    ! it here so the namelist reads straight into the single source of truth.
    use profile_data_mod, only: boole_kim_reff_coords

    implicit none

    ! Module-level storage for delta-f weight-evolution inputs.
    ! These are read by read_rmp_response_currents_inp_into_type and
    ! consumed by calc_rmp_response_currents after grid build.
    character(len=512), public :: profile_dir          = './profiles'
    character(len=512), public :: equil_mapping_file   = './flux_functions.dat'
    logical,            public :: boole_constant_delta_B_r = .true.
    real(dp),           public :: delta_B_r_const    = 0.0_dp
    integer,            public :: pert_m_fourier       = 0
    integer,            public :: pert_n_fourier       = 0
    character(len=512), public :: delta_B_r_file     = ''
    ! Step-function Br mode: Br = delta_B_r_const in [centre-hw, centre+hw], 0 outside.
    ! step_center_reff and step_halfwidth_reff are in r_eff [cm].
    ! r_eff is the toroidal-flux-based effective radius (NOT the geometric minor
    ! radius r_geom).  Read r_eff from flux_functions.dat col1 (0-indexed).
    logical,  public :: boole_step_delta_B_r   = .false.
    real(dp), public :: step_center_reff       = 0.0_dp  ! r_eff [cm] — NOT r_geom
    real(dp), public :: step_halfwidth_reff    = 0.0_dp  ! r_eff [cm] — NOT r_geom
    ! Perturbed electrostatic E_perp from KIM (E_perp.dat).
    ! boole_e_perp enables the c*dE^s/B0 term in wdot (Albert 2016 Eq.4,
    ! extended to include E×B radial drift from the RMP electrostatic response).
    ! e_perp_file must be KIM's E_perp.dat (4-column: r_eff[cm], Re, Im, |mag|).
    ! Must be used together with one of the Br perturbation modes (shares m,n).
    ! r_eff here is the toroidal-flux-based effective radius (NOT r_geom).
    logical,            public :: boole_e_perp  = .false.
    character(len=512), public :: e_perp_file   = ''   ! path to KIM E_perp.dat
    integer,            public :: species_for_delta_f  = 1
    ! Diagnostic: if .true., skip the exp(i*(m theta + n phi)) factor in
    ! the constant-amplitude perturbation and return delta_B_r_const as
    ! a truly uniform (axisymmetric) real-valued perturbation.
    logical,            public :: boole_skip_phase_for_test = .false.

    ! Regularisation parameters (Albert 2016 Eq. 4 with linear damping).
    ! The damping rate is nu_r = nu_r_frac * nu_c per marker (still tied
    ! to the local collision frequency, since nu_r is a physics quantity),
    ! and the damping switches on at M * tau_c. The total trace time is
    ! now controlled by the global namelist parameter time_step -- every
    ! marker traces for the same wall-clock interval.
    real(dp),           public :: nu_r_frac                = 0.5_dp
    integer,            public :: m_collision_times_reg_on = 3
    ! Coulomb logarithm for the local collision frequency.
    real(dp),           public :: coulomb_log              = 17.0_dp
    ! Multiplicative scaling on the collisional rates (dpp, dhh, fpeff) used
    ! inside stost. Used to benchmark against external collision-frequency
    ! conventions (e.g. KIM nu_e). Default 1.0 = no change.
    real(dp),           public :: nu_scale_factor          = 1.0_dp
    ! Use an external collision frequency table for the OU operator
    ! (iswmode=5). When enabled, stost gets nu_step from the splined
    ! nu_e(s) loaded from kim_nu_file instead of 4*dhh.  FP runs are
    ! unaffected (the override only triggers in the iswmode=5 branch).
    logical,            public :: boole_use_kim_nu         = .false.
    character(len=256), public :: kim_nu_file              = ''

    ! Radially-varying anomalous diffusion D_a(r_eff) profile.
    ! When da_profile_file is non-empty, the profile is loaded and the
    ! per-marker D_a at each step is evaluated from the spline instead of
    ! using the scalar in%anomalous_diffusion_coefficient.
    ! da_scale_factor multiplies the profile values (use for sensitivity scans).
    ! To enable the profile kick, set anomalous_diffusion_coefficient > 0
    ! in the namelist (acts as a master enable; its exact value is ignored
    ! once da_profile_file is non-empty) OR the profile gate supersedes it.
    character(len=256), public :: da_profile_file          = ''
    real(dp),           public :: da_scale_factor          = 1.0_dp

    ! Diagnostic: when true, bypass the delta-f weight evolution entirely
    ! and hold w(t) = 1 (real, constant) for the whole trace. Useful for
    ! tracer-style runs where each marker should deposit its raw v_par
    ! contribution at unit weight, decoupled from the Albert delta-f
    ! amplitude/regularisation machinery. Set nu_r_frac = 0 alongside.
    logical,            public :: boole_constant_unit_weight = .false.

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

    ! Maximum number of respawn attempts per particle. When > 0 and delta-f
    ! is on, a particle that exits the computational domain (ind_tetr=-1)
    ! has its starting position re-drawn (within the bias window intersected
    ! with s_outer_cut) and continues being traced; deposits accumulate
    ! across attempts and the per-particle time average uses the SUM of
    ! actual t_confined across attempts. Default = 0 (no respawn).
    integer,            public :: n_respawn_max = 0

    ! When true, replace the rejection-sampled bias of starting positions
    ! with one marker per radial s-layer (theta and phi uniform in their
    ! valid ranges). The number of placed markers is the count of layer
    ! midpoints falling in [s_inner_sample, s_outer_sample] -- so set
    ! n_particles in the input to at least this count (typically ~ n1).
    ! Layers in excess of n_particles slots are skipped; unused slots are
    ! marked lost.
    logical,            public :: boole_equidistant_s_sampling = .false.

    ! Within each radial layer of the equidistant-in-s spawner, place
    ! markers on a stratified theta grid: marker k of K in a layer draws
    ! theta uniformly inside the bin [(k-1)*2pi/K, k*2pi/K). Reduces
    ! poloidal-phase clustering at small K. Only takes effect when
    ! boole_equidistant_s_sampling = .true.
    logical,            public :: boole_stratify_theta = .false.

    ! Same idea for phi. When BOTH stratification flags are on, the phi
    ! bin assignment for marker k uses a per-layer random permutation
    ! (Fisher-Yates) of {1..K} -- a Latin-hypercube layout so that the
    ! (theta, phi) pairs aren't locked to a diagonal.
    logical,            public :: boole_stratify_phi = .false.

    ! Diagnostic: when true, dump marker n=1's trajectory (R, phi, Z,
    ! vpar, ind_tetr) at sparse intervals to traj_n1.dat for comparison
    ! of orbit-derived q vs qsaf from flux_functions.dat. Negligible
    ! overhead since only one marker is logged.
    logical,            public :: boole_dump_orbit_n1 = .false.
    integer,            public :: orbit_dump_stride  = 50
    integer,            public :: traj_dump_unit     = 0
    ! Persistent push-step counter; used by the trajectory dump to stride.
    integer,            public :: traj_step_count    = 0

    ! Collision-event diagnostic for marker n=1: writes one row per call to
    ! carry_out_collisions to coll_n1.dat, with the pre-push tetra, the
    ! step size dt to the next collision, |v|, etc. Counters (event total
    ! and dt sum) accumulate on every event so end-of-run statistics are
    ! exact regardless of stride. The stride only governs which rows hit
    ! disk -- coll_n1.dat gets every coll_dump_stride-th event.
    logical, public :: boole_dump_collisions_n1 = .false.
    integer, public :: coll_dump_unit          = 0
    integer, public :: coll_dump_stride        = 1
    integer, public :: coll_event_count        = 0
    real(dp), public :: coll_dt_sum            = 0.0_dp
    real(dp), public :: coll_dist_sum          = 0.0_dp

    ! Point-source spawn position when boole_point_source = .true.,
    ! interpreted in the integrator's native (a, b, c) chart. For
    ! grid_kind=2 / cylindrical that's (R, phi, Z) in cm/rad/cm. If left
    ! at its sentinel default (all zeros) the legacy fallback in
    ! set_starting_positions_rmp is used.
    real(dp), public :: point_source_x(3) = 0.0_dp

    ! Override marker n=1's pitch (v_par/v) to a fixed value. Useful for
    ! single-marker diagnostics where you want a deterministic orbit
    ! (pitch=1 -> pure field-line follower, v_perp=0). When .false. the
    ! pitch sampled from the normal distribution is kept.
    logical,  public :: boole_force_marker1_pitch = .false.
    real(dp), public :: marker1_pitch_value       = 1.0_dp

    ! Population filter: keep only passing or only trapped markers.
    ! 'none' (default), 'passing', or 'trapped'. Filtered-out markers
    ! get start%lost = .true. and are skipped during tracing. Useful for
    ! decomposing the screening current into its passing vs trapped
    ! contributions (Fortran sums them; their plotted magnitudes add up
    ! to the unfiltered result with the same n_particles divisor).
    character(len=16),  public :: trapping_filter_mode = 'none'

    ! Initial-energy distribution for the markers.
    !   'mono'        : delta function at in%energy_eV
    !   'uniform'     : flat in [0, e_cutoff_factor * in%energy_eV]
    !   'maxwellian'  : p(E) ∝ sqrt(E) * exp(-E/T) on the same range,
    !                   T = in%energy_eV (velocity-Jacobian-weighted Maxwellian)
    character(len=32), public :: energy_dist_kind = 'maxwellian'
    real(dp),          public :: e_cutoff_factor  = 5.0_dp
    ! Module-private temperature plumbed into pdf_maxwellian_energy.
    real(dp)                  :: T_sample_eV     = 0.0_dp

    ! Per-particle regularisation storage. Allocated alongside weights%w
    ! when boole_delta_f is on. tau_c is the local collision time at the
    ! starting position; t_reg_on the switch-on time of the damping;
    ! nu_r the damping rate.
    real(dp), allocatable, public :: tau_c(:,:), t_reg_on(:,:)
    real(dp), allocatable, public :: nu_r(:,:)

contains

subroutine read_rmp_response_currents_inp_into_type

    use gorilla_applets_types_mod, only: in
    use collis_ions,               only: collis_nu_scale_factor => nu_scale_factor

    real(dp) :: time_step, energy_eV, n_particles, density
    real(dp) :: anomalous_diffusion_coefficient
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, &
               boole_write_exit_data, boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, &
               boole_eliminate_particles_outside_flux, boole_delta_f
    integer :: i_integrator_type, seed_option, n_species
    real(dp) :: flux_threshold_for_elimination
    ! Collision-operator mode for carry_out_collisions:
    !   1 = full Coulomb (stost iswmode=1, current default)
    !   5 = Euler-Maruyama Ornstein-Uhlenbeck step on v_par (stost iswmode=5)
    ! Default 1 keeps legacy runs byte-identical.
    integer :: i_collision_mode = 1

    integer :: s_inp_unit

    ! Note: energy_dist_kind, e_cutoff_factor are module-level publics with
    ! defaults at declaration; the NAMELIST references them directly so the
    ! namelist key matches the module variable name 1:1.
    NAMELIST /rmp_response_currents_nml/ time_step, energy_eV, n_particles, boole_squared_moments, boole_point_source, &
    & boole_collisions, boole_precalc_collisions, density, boole_refined_sqrt_g, &
    & energy_dist_kind, e_cutoff_factor, &
    & boole_linear_density_simulation, boole_antithetic_variate, boole_linear_temperature_simulation, i_integrator_type, &
    & seed_option, boole_write_vertex_indices, boole_write_vertex_coordinates, boole_write_prism_volumes, &
    & boole_write_refined_prism_volumes, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_species, &
    & boole_eliminate_particles_outside_flux, flux_threshold_for_elimination, boole_delta_f, &
    & profile_dir, equil_mapping_file, boole_constant_delta_B_r, delta_B_r_const, &
    & pert_m_fourier, pert_n_fourier, delta_B_r_file, &
    & boole_step_delta_B_r, step_center_reff, step_halfwidth_reff, species_for_delta_f, &
    & nu_r_frac, m_collision_times_reg_on, coulomb_log, nu_scale_factor, &
    & boole_use_kim_nu, kim_nu_file, &
    & boole_kim_reff_coords, &
    & da_profile_file, da_scale_factor, &
    & boole_constant_unit_weight, &
    & boole_compute_n_modes_dft, s_outer_cut, boole_skip_phase_for_test, &
    & boole_e_perp, e_perp_file, &
    & s_inner_sample, s_outer_sample, n_respawn_max, boole_equidistant_s_sampling, &
    & boole_stratify_theta, boole_stratify_phi, &
    & boole_dump_orbit_n1, orbit_dump_stride, trapping_filter_mode, &
    & point_source_x, boole_force_marker1_pitch, marker1_pitch_value, &
    & boole_dump_collisions_n1, coll_dump_stride, i_collision_mode, &
    & anomalous_diffusion_coefficient

    ! Default: no anomalous transport (D_anom = 0 disables the kick).
    anomalous_diffusion_coefficient = 0.0_dp

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
    in%i_collision_mode = i_collision_mode
    in%boole_refined_sqrt_g = boole_refined_sqrt_g
    ! Sanity-fix e_cutoff_factor if user passed a non-positive value.
    if (e_cutoff_factor <= 0.0_dp) e_cutoff_factor = 5.0_dp
    ! Note: in%boole_monoenergetic is no longer used in the rmp applet --
    ! the new energy_dist_kind module variable takes its place.
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
    in%anomalous_diffusion_coefficient = anomalous_diffusion_coefficient

    ! Propagate nu_scale_factor to collis_ions module variable used inside stost.
    collis_nu_scale_factor = nu_scale_factor

    print '(a, f10.4)', ' nu_scale_factor (stost rate multiplier) = ', collis_nu_scale_factor

end subroutine read_rmp_response_currents_inp_into_type

! ====================================================================
! Allocates the per-particle regularisation arrays. Initialised to 0.
! ====================================================================
subroutine allocate_delta_f_per_particle(n_particles, n_species)

    integer, intent(in) :: n_particles, n_species

    if (.not. allocated(tau_c))      allocate(tau_c(n_particles, n_species))
    if (.not. allocated(t_reg_on))   allocate(t_reg_on(n_particles, n_species))
    if (.not. allocated(nu_r))       allocate(nu_r(n_particles, n_species))

    tau_c      = 0.0_dp
    t_reg_on   = 0.0_dp
    nu_r       = 0.0_dp

end subroutine allocate_delta_f_per_particle

! ====================================================================
! Initial-condition setup, copied from utils_data_pre_and_post_processing_mod
! so the rmp applet owns the chain and can iterate freely (in particular on
! the energy distribution).
! ====================================================================
subroutine calc_starting_conditions_rmp_response_currents(verts)

    use gorilla_applets_types_mod, only: in, start

    real(dp), dimension(:,:), allocatable, intent(out), optional :: verts

    call set_verts_and_coordinate_limits_rmp(verts)

    call allocate_start_type_rmp
    call allocate_weights_rmp
    call set_starting_positions_rmp()
    call set_rest_of_start_type_rmp()

end subroutine calc_starting_conditions_rmp_response_currents

! --------------------------------------------------------------------
subroutine set_verts_and_coordinate_limits_rmp(verts)

    use tetra_physics_mod, only: coord_system, tetra_physics
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, nvert
    use tetra_grid_settings_mod, only: grid_size
    use magdata_in_symfluxcoor_mod, only : raxis, zaxis
    use gorilla_applets_types_mod, only: g

    real(dp), dimension(:,:), allocatable, intent(out), optional :: verts
    real(dp), dimension(:,:), allocatable :: verts_local
    integer :: i

    g%ind_a = 1 !(R in cylindrical coordinates, s in flux coordinates)
    g%ind_b = 2 !(phi in cylindrical and flux coordinates)
    g%ind_c = 3 !(z in cylindrical coordinates, theta in flux coordinates)
    if (coord_system.eq.2) then
        g%ind_b = 3
        g%ind_c = 2
    endif

    allocate(verts_local(3,nvert))
    if (coord_system.eq.1) verts_local = verts_rphiz
    if (coord_system.eq.2) verts_local = verts_sthetaphi

    g%amin = minval(verts_local(g%ind_a,:))
    g%amax = maxval(verts_local(g%ind_a,:))
    g%cmin = minval(verts_local(g%ind_c,:))
    g%cmax = maxval(verts_local(g%ind_c,:))

    g%raxis = raxis
    g%zaxis = zaxis

    g%dist_from_o_point_within_grid = 0.0_dp
    do i = 1,3*grid_size(3)
        g%dist_from_o_point_within_grid = max(g%dist_from_o_point_within_grid, &
                                              1.1_dp*sqrt((tetra_physics(i)%x1(1)-raxis)**2 + (tetra_physics(i)%x1(3)-zaxis)**2))
    enddo

    if (present(verts)) call move_alloc(verts_local, verts)

end subroutine set_verts_and_coordinate_limits_rmp

! --------------------------------------------------------------------
subroutine allocate_start_type_rmp

    use gorilla_applets_types_mod, only: start, in
    use gorilla_applets_settings_mod, only: i_option

    allocate(start%x(3,in%num_particles,in%n_species))
    allocate(start%pitch(in%num_particles,in%n_species))
    allocate(start%energy(in%num_particles,in%n_species))
    allocate(start%jperp(in%num_particles,in%n_species))
    allocate(start%lost(in%num_particles,in%n_species))
    allocate(start%particle_charge(in%n_species))
    allocate(start%particle_mass(in%n_species))
    allocate(start%cm_over_e(in%n_species))
    allocate(start%t(in%n_species))
    allocate(start%v0(in%n_species))

end subroutine allocate_start_type_rmp

! --------------------------------------------------------------------
subroutine allocate_weights_rmp

    use gorilla_applets_types_mod, only: in, weights

    allocate(weights%w(in%num_particles,in%n_species))
    weights%w = (0.0_dp, 0.0_dp)
    if (in%boole_delta_f) then
        allocate(weights%original(in%num_particles,in%n_species))
        weights%original = (0.0_dp, 0.0_dp)
    end if

end subroutine allocate_weights_rmp

! --------------------------------------------------------------------
subroutine set_starting_positions_rmp()

    use gorilla_applets_types_mod, only: in, start, g
    use tetra_physics_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind
    use constants, only: pi
    use marker_distribution_mod, only: pdf_flat, init_distribution_3d, sample_array_3d

    real(dp) :: xmin(3), xmax(3)
    integer :: i_species

    !compute starting conditions
    if (in%boole_point_source) then
        if (any(point_source_x /= 0.0_dp)) then
            xmin = point_source_x
        elseif (grid_kind.eq.2) then
            xmin = [209.0_dp, 0.01_dp, 10.0_dp]
        elseif (grid_kind.eq.4) then
            xmin = [205.0_dp, 0.0_dp, 0.0_dp]
        endif
        xmax = xmin
        if (coord_system.eq.2) print*, 'error: point source is only implemented for cylindrical coordinate system'
    else
        ! Set bounds for uniform sampling
        xmin(g%ind_a) = g%amin
        xmax(g%ind_a) = g%amax
        xmin(g%ind_b) = 0.0_dp
        xmax(g%ind_b) = 2*pi
        xmin(g%ind_c) = g%cmin
        xmax(g%ind_c) = g%cmax
    endif

    ! Initialize position distribution if not already done
    if (.not. start%dist_position%initialized) then
        call init_distribution_3d(start%dist_position, pdf_flat, xmin, xmax)
    endif

    ! Sample positions for all species
    do i_species = 1, in%n_species
        call sample_array_3d(start%dist_position, start%x(:,:,i_species))
    enddo

end subroutine set_starting_positions_rmp

! --------------------------------------------------------------------
subroutine set_rest_of_start_type_rmp()

    use gorilla_applets_types_mod, only: in, start
    use tetra_physics_mod, only: cm_over_e, particle_charge, particle_mass
    use gorilla_applets_settings_mod, only: i_option
    use constants, only: echarge, ame, clight
    use constants, only: ev2erg
    use marker_distribution_mod, only: pdf_flat, init_distribution_1d, sample_array_1d
    use profile_data_mod, only: load_da_profile

    integer :: i_species

    ! Initialize lambda (pitch angle) distribution if not already done
    if (.not. start%dist_lambda%initialized) then
        call init_distribution_1d(start%dist_lambda, pdf_flat, -1.0_dp, 1.0_dp)
    endif

    ! Initialize energy distribution if not already done.
    ! Dispatch on energy_dist_kind:
    !   'mono'        : delta function at in%energy_eV
    !   'uniform'     : flat in [0, e_cutoff_factor * in%energy_eV]
    !   'maxwellian'  : sqrt(E)*exp(-E/T) on [0, e_cutoff_factor * T], T = in%energy_eV
    if (.not. start%dist_energy%initialized) then
        select case(trim(energy_dist_kind))
        case('mono')
            call init_distribution_1d(start%dist_energy, pdf_flat, &
                                      in%energy_eV, in%energy_eV)
        case('uniform')
            call init_distribution_1d(start%dist_energy, pdf_flat, &
                                      0.0_dp, e_cutoff_factor * in%energy_eV)
        case('maxwellian')
            T_sample_eV = in%energy_eV       ! plumb T into the pdf function
            call init_distribution_1d(start%dist_energy, pdf_maxwellian_energy, &
                                      0.0_dp, e_cutoff_factor * in%energy_eV)
        case default
            print *, "Error: unknown energy_dist_kind = '", trim(energy_dist_kind), &
                     "' (expected 'mono', 'uniform', or 'maxwellian')"
            error stop
        end select
    endif

    ! Sample pitch and energy for all species
    do i_species = 1, in%n_species
        call sample_array_1d(start%dist_lambda, in%num_particles, start%pitch(:,i_species))
        call sample_array_1d(start%dist_energy, in%num_particles, start%energy(:,i_species))
    enddo

    if (in%boole_antithetic_variate) then
        start%x(:,1:in%num_particles:2,:) = start%x(:,2:in%num_particles:2,:)
        start%pitch(1:in%num_particles:2,:) = -start%pitch(2:in%num_particles:2,:)
        start%energy(1:in%num_particles:2,:) = start%energy(2:in%num_particles:2,:)
    endif

    ! Diagnostic override of marker n=1's pitch. Controlled by the
    ! namelist flag boole_force_marker1_pitch and the value
    ! marker1_pitch_value. Default OFF -- the normally sampled pitch is
    ! kept. Useful for single-marker tests where a deterministic orbit
    ! is wanted (e.g. pitch=1 for a pure field-line follower).
    if (boole_force_marker1_pitch .and. in%num_particles >= 1) then
        do i_species = 1, in%n_species
            start%pitch(1, i_species) = marker1_pitch_value
        end do
    end if

    start%particle_charge = particle_charge
    start%particle_mass = particle_mass
    start%cm_over_e = cm_over_e
    start%t = in%time_step
    if (i_option.eq.12) then
        start%particle_charge(2) = -echarge
        start%particle_mass(2) = ame
        start%cm_over_e(2) = -clight*ame/echarge
        start%t(2) = in%time_step/42.0_dp
    endif

    start%v0 = sqrt(2.0_dp*in%energy_eV*ev2erg/start%particle_mass)
    start%lost = .false.

    ! Debug: dump sampled initial energies (one per line, eV) for sampler verification.
    block
        integer :: u, ip
        open(newunit=u, file='initial_energies.dat', status='replace', action='write')
        do ip = 1, in%num_particles
            write(u,'(es16.8)') start%energy(ip, 1)
        enddo
        close(u)
    end block

    ! Load radially-varying D_a(r_eff) profile if requested.
    if (len_trim(da_profile_file) > 0) then
        call load_da_profile(trim(da_profile_file))
    end if

end subroutine set_rest_of_start_type_rmp

! ====================================================================
subroutine parallelised_particle_pushing_rmp_response_currents(species, n_particles_in)

    use gorilla_applets_types_mod, only: counter, c, in, time_t, moment_specs, counter_t, particle_status_t, start, s
    use tetra_grid_mod, only: ntetr
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use utils_parallelised_particle_pushing_mod, only: print_progress, handle_lost_particles, add_local_tetr_moments_to_output, &
        add_local_counter_to_counter, initialise_loop_variables, carry_out_collisions, update_exit_data, update_start_type, &
        initialise_seed_for_random_numbers_for_each_thread
    use find_tetra_mod, only: find_tetra
    ! Anomalous-transport kick: re-use the existing applet's perpendicular
    ! random-walk displacement. Enabled when in%anomalous_diffusion_coefficient > 0
    ! or da_profile_file is loaded.
    use anomalous_transport_displacement_mod, only: anomalous_transport_displacement
    use profile_data_mod, only: da_profile_loaded, eval_da_profile

    integer, intent(in)                               :: species
    integer, intent(in), optional                     :: n_particles_in
    integer                                           :: kpart, iantithetic, ind_tetr, iface, n_particles
    integer                                           :: p, l, n, i, i_total, n_respawn_used, n_respawn_total
    integer                                           :: n_truly_lost
    real(dp), dimension(3)                            :: x
    real(dp)                                          :: vpar, vperp, t_tot, trace_time_n, t_actual_n
    type(time_t)                                      :: t
    type(counter_t)                                   :: local_counter
    type(particle_status_t)                           :: particle_status
    complex(dp), dimension(:,:), allocatable          :: local_tetr_moments
    ! Per-particle workspace for the delta-f time average. Each marker's
    ! deposits accumulate here across ALL respawn attempts (we do NOT reset
    ! the workspace on respawn). At the end of all attempts we divide by
    ! t_actual_n -- the SUM of t_confined across attempts -- and fold into
    ! local_tetr_moments. So even a marker that is killed and never
    ! respawned contributes its short-trace deposit divided by its actual
    ! short trace time, which is statistically equivalent to its share of a
    ! full trace (rather than being suppressed by the full intended T_n).
    complex(dp), dimension(:,:), allocatable          :: particle_tetr_moments
    real(dp)                                          :: da_local
    integer                                           :: n_da_sub, i_da_sub
    logical                                           :: thread_flag = .true.
    logical                                           :: respawn_success

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

    n_respawn_total = 0
    n_truly_lost    = 0

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(counter, kpart, species, in, c, iantithetic, start, s, n_particles, moment_specs, ntetr, n_respawn_max, &
    !$OMP&        boole_dump_orbit_n1, traj_dump_unit, orbit_dump_stride, traj_step_count, &
    !$OMP&        boole_dump_collisions_n1, coll_dump_unit, coll_dump_stride, &
    !$OMP&        coll_event_count, coll_dt_sum, coll_dist_sum, &
    !$OMP&        da_profile_loaded, da_scale_factor) &
    !$OMP& REDUCTION(+:t_tot, n_respawn_total, n_truly_lost) &
    !$OMP& PRIVATE(p, l, n, i, i_total, n_respawn_used, x, vpar, vperp, t, ind_tetr, iface, local_tetr_moments, local_counter, particle_status, trace_time_n, particle_tetr_moments, t_actual_n, respawn_success, da_local, n_da_sub, i_da_sub) &
    !$OMP& FIRSTPRIVATE(thread_flag)

    if (omp_get_thread_num().eq.0) print*, 'Number of threads: ', omp_get_num_threads()

    ! Per-thread allocation of the per-particle workspace. Allocatables in
    ! PRIVATE start unallocated on each thread, so this must happen inside
    ! the parallel region. Only used in the delta-f path; harmless to
    ! allocate either way.
    if (in%boole_delta_f) then
        allocate(particle_tetr_moments(moment_specs%n_moments, ntetr))
    end if

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

            ! Per-particle accumulators that persist across respawn attempts.
            if (in%boole_delta_f) particle_tetr_moments = (0.0_dp, 0.0_dp)
            n_respawn_used = 0

            ! One-shot init: place the marker, zero t%confined, zero the
            ! delta-f weight (via particle_status%initialized = .false.
            ! triggering the init block in orbit_timestep_rmp_response_currents).
            call initialise_loop_variables(l, n, local_counter, particle_status, t, local_tetr_moments, x, vpar, vperp, species)

            ! Total trace time is the global time_step for every marker
            ! (single source of truth). The per-marker damping rate nu_r
            ! is still local to each marker's spawn s; it gets
            ! initialised lazily inside the first orbit_timestep call,
            ! and is REFRESHED at every respawn site below.
            trace_time_n = start%t(species)

            i = 0

            do while (t%confined.lt.trace_time_n)
                i = i + 1

                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, t, x, vpar, vperp, ind_tetr, iface, species, iswmode_in=in%i_collision_mode)
                    t%step = t%step / start%v0(species)
                    ! Collision-event diagnostics for marker n=1:
                    ! accumulate counters on every event (so end-of-run
                    ! stats are exact), but write coll_n1.dat only every
                    ! coll_dump_stride-th event to keep file size sane.
                    if (boole_dump_collisions_n1 .and. n == 1 .and. ind_tetr /= -1) then
                        !$omp critical (coll_dump)
                        coll_event_count = coll_event_count + 1
                        coll_dt_sum      = coll_dt_sum   + t%step
                        coll_dist_sum    = coll_dist_sum + sqrt(vpar**2 + vperp**2) * t%step
                        if (coll_dump_unit /= 0 .and. mod(coll_event_count, max(coll_dump_stride,1)) == 0) then
                            write(coll_dump_unit, '(i12, es16.8, i10, 3es16.8, 3es16.8, es16.8)') &
                                coll_event_count, t%confined, ind_tetr, &
                                x(1), x(2), x(3), vpar, vperp, &
                                sqrt(vpar**2 + vperp**2), t%step
                        end if
                        !$omp end critical (coll_dump)
                    end if
                else
                    t%step = trace_time_n - t%confined
                endif

                if (in%boole_delta_f) then
                    call orbit_timestep_rmp_response_currents(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                                              particle_tetr_moments, local_counter, species, trace_time_n)
                else
                    call orbit_timestep_rmp_response_currents(x, vpar, vperp, t, particle_status, ind_tetr, iface, n, &
                                                              local_tetr_moments, local_counter, species, trace_time_n)
                end if

                t%confined = t%confined + t%step - t%remain
                t_tot = t_tot + t%step - t%remain

                ! Anomalous-transport kick (perpendicular random walk + correction
                ! velocity), re-used from anomalous_transport_displacement_mod.
                ! Disabled when neither the scalar nor a profile is active. The weight
                ! weights%w(n, species) is intentionally NOT reset here -- the
                ! marker carries its coherent weight to the displaced location,
                ! where wdot_s at the new (s, theta, phi) drives the subsequent
                ! evolution. This matches the soft-respawn convention used for
                ! out-of-domain re-entries.
                if ((in%anomalous_diffusion_coefficient > 0.0_dp .or. da_profile_loaded) &
                    .and. ind_tetr /= -1) then
                    if (da_profile_loaded) then
                        da_local = eval_da_profile(eval_s_local(ind_tetr, x)) * da_scale_factor
                    else
                        da_local = in%anomalous_diffusion_coefficient
                    end if
                    if (da_local > 0.0_dp) then
                        ! Sub-step when a single orbit step would displace > epsilon=0.1 cm.
                        ! n_da_sub is chosen so each sub-step displacement <= epsilon.
                        ! t%step is partitioned evenly so total diffusion per orbit step
                        ! = sqrt(2*D*t%step) regardless of n_da_sub.
                        n_da_sub = max(1, nint(2.0_dp * da_local * t%step / 1.0d-2))
                        t%step_anomalous_transport = t%step / real(n_da_sub, dp)
                        do i_da_sub = 1, n_da_sub
                            call anomalous_transport_displacement(x, ind_tetr, iface, &
                                t%step_anomalous_transport, vpar, vperp, da_local)
                            if (ind_tetr == -1) exit
                        end do
                    end if
                end if

                if (ind_tetr == -1) then
                    ! Soft respawn: redraw only the spatial position. Keep
                    ! t%confined, weights%w, vpar, vperp, and the running
                    ! deposit. Refresh nu_r/tau_c at the new local s.
                    ! Tried hard reset (w=0 on respawn) -- concentrates the
                    ! (0,0) Fourier residual at the resonance shell, ~3x
                    ! worse than soft. Soft is the only sensible option.
                    if (.not. in%boole_delta_f .or. n_respawn_used >= n_respawn_max) exit
                    call respawn_particle_in_window(n, species, respawn_success)
                    if (.not. respawn_success) exit
                    x = start%x(:, n, species)
                    call find_tetra(x, vpar, vperp, ind_tetr, iface)
                    if (ind_tetr == -1) exit
                    call init_regularisation_for_particle(n, ind_tetr, x, vpar, vperp, species)
                    particle_status%lost = .false.
                    particle_status%exit = .false.
                    ! particle_status%initialized stays .true. so the
                    ! delta-f weight is NOT zeroed here.
                    n_respawn_used  = n_respawn_used  + 1
                    n_respawn_total = n_respawn_total + 1
                end if
            enddo

            i_total    = i
            t_actual_n = t%confined

            if (ind_tetr == -1 .and. t%confined < trace_time_n) then
                call handle_lost_particles(local_counter, particle_status%lost)
                n_truly_lost = n_truly_lost + 1
            end if

            ! Per-particle time average using the SUM of t_confined over all
            ! attempts (= the actual integration time we have data for).
            ! This makes killed-and-not-respawned markers contribute their
            ! short trace at full weight, instead of being suppressed by
            ! division by the full intended T_n they never reached.
            if (in%boole_delta_f .and. t_actual_n > 0.0_dp) then
                local_tetr_moments = local_tetr_moments &
                                   + particle_tetr_moments / t_actual_n
            end if

            !$omp critical
            counter%integration_steps = counter%integration_steps + i_total
            c%maxcol = max(dble(i_total)/dble(c%randcoli), c%maxcol)
            call add_local_counter_to_counter(local_counter)
            !$omp end critical

            call update_exit_data(particle_status%lost, t_actual_n, x, vpar, vperp, i_total, n, species_in=species, ind_tetr=ind_tetr)
            call update_start_type(x, vpar, vperp, n, species, ind_tetr)
        enddo

        !$omp critical
        call add_local_tetr_moments_to_output(local_tetr_moments, species)
        !$omp end critical
    enddo
    !$OMP END DO

    if (allocated(particle_tetr_moments)) deallocate(particle_tetr_moments)
    !$OMP END PARALLEL

    print*, 'Total tracing time / number of particles: ', t_tot/n_particles, 's'
    if (in%boole_delta_f .and. n_respawn_max > 0) then
        print*, 'Respawn attempts used (across all particles): ', n_respawn_total
        print*, 'Particles still lost after respawn budget   : ', n_truly_lost
    end if

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

    real(dp), dimension(3)                       :: z_save, x_new, x_pre_push
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
            ! Per-particle regularisation parameters from local profile values.
            ! In unit-weight mode we still need trace_time to be set so the
            ! trace duration matches the standard delta-f run; nu_r and
            ! t_reg_on go unused because update_delta_f_weight is a no-op.
            call init_regularisation_for_particle(n, ind_tetr, x, vpar, vperp, species)
            if (boole_constant_unit_weight) then
                ! Tracer mode: hold w(t) = 1 + 0i for the whole trace.
                weights%w(n, species) = (1.0_dp, 0.0_dp)
                weights%original(n, species) = (1.0_dp, 0.0_dp)
            else
                ! Always start from w(0) = 0 (paper convention).
                weights%w(n, species) = (0.0_dp, 0.0_dp)
                weights%original(n, species) = (0.0_dp, 0.0_dp)
            end if
        end if
        particle_status%initialized = .true.
    endif

    if (t%step.eq.0.0_dp) return
    if (particle_status%initialized) z_save = x - tetra_physics(ind_tetr)%x1
    call initialize_constants_of_motion(vperp, z_save, ind_tetr, perpinv)

    t%remain = t%step
    boole_t_finished = .false.
    local_counter%tetr_pushings = local_counter%tetr_pushings - 1

    ! Safe defaults so axis-recovery can read these on the very first iteration
    ! even if ind_tetr is already -1 when orbit_timestep is entered.
    x_pre_push = x
    ind_tetr_save = ind_tetr

    do
        local_counter%tetr_pushings = local_counter%tetr_pushings + 1

        if (ind_tetr.eq.-1) then
            if ((grid_kind.eq.2).or.(grid_kind.eq.3)) then
                call identify_particles_entering_annulus(x, local_counter, boole_lost_inside)
                if (boole_lost_inside) then
                    x_new = 3*(/g%raxis, x(2), g%zaxis/) - 2*x
                    vperp = vperp_func(z_save, perpinv, ind_tetr_save)
                    call find_tetra(x_new, vpar, vperp, ind_tetr, iface)
                    if (ind_tetr.ne.-1) then
                        x = x_new
                    else
                        ! Axis reflection failed (reflected position outside mesh).
                        ! Exit with ind_tetr=-1 so the outer loop can respawn the marker.
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
        x_pre_push = x

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
        ! Midpoint quadrature in cell ind_tetr_save: wdot_s is evaluated
        ! at the spatial midpoint of (x_pre_push, x_post_push), using the
        ! tetra_physics of the cell the marker just traversed. Previously
        ! we evaluated at x_post_push with the NEW cell's physics — a
        ! left-Riemann sample at the wrong boundary.
        if (in%boole_delta_f .and. ind_tetr_save /= -1) then
            call update_delta_f_weight(n, ind_tetr_save, &
                                       0.5_dp * (x_pre_push + x), &
                                       vpar, vperp, t_pass, t, species)
        endif

        call update_local_tetr_moments(local_tetr_moments, ind_tetr_save, n, optional_quantities, species)
        if ((grid_kind.eq.2).or.(grid_kind.eq.3)) call compute_radial_fluxes(ind_tetr_save, ind_tetr, x)

        ! Diagnostic: dump marker n=1 trajectory for orbit-q comparison.
        ! Use ind_tetr (the tetra containing the post-push x), NOT
        ! ind_tetr_save (the pre-push tetra) -- otherwise barycentric
        ! interpolation in the post-processor extrapolates outside the
        ! tetra and gives spurious s excursions.
        if (boole_dump_orbit_n1 .and. n == 1 .and. traj_dump_unit /= 0 .and. ind_tetr /= -1) then
            traj_step_count = traj_step_count + 1
            if (mod(traj_step_count, orbit_dump_stride) == 0) then
                !$omp critical (traj_dump)
                write(traj_dump_unit, '(es16.8, 3es16.8, es16.8, 4es16.8, i10, i10)') &
                    t%confined + t%step - t%remain, x(1), x(2), x(3), vpar, &
                    real(weights%w(n, species), dp), aimag(weights%w(n, species)), &
                    real(weights%w(n, species), dp) * vpar, &
                    aimag(weights%w(n, species)) * vpar, &
                    ind_tetr, iper_phi
                !$omp end critical (traj_dump)
            end if
        end if

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
    use tetra_physics_mod, only: tetra_physics, metric_determinant
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
    else if (coord_system == 2) then
        ! Symmetry-flux coordinates: the spatial Jacobian is sqrt(g), NOT the
        ! cylindrical 'R = r + x1(1)' formula below. Use the analytic metric
        ! determinant (clean and finite for grid_kind=5; the refined sqrt_g
        ! formula above can hit a curlA-based NaN in degenerate seam tetras).
        J_x = metric_determinant(ind_tetr, z_save + tetra_physics(ind_tetr)%x1)
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

    ! Tracer mode: weight is fixed at 1, no evolution.
    if (boole_constant_unit_weight) return

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
    use constants, only: ev2erg, clight
    use profile_data_mod, only: eval_profiles, profile_values_t
    use perturbation_field_mod, only: eval_delta_B_s, eval_delta_E_s

    integer,     intent(in)  :: ind_tetr, species
    real(dp),    intent(in)  :: x(3), vpar, vperp
    complex(dp), intent(out) :: wdot_s

    real(dp)    :: z_cell(3), s_loc, theta_loc, phi_loc
    real(dp)    :: B0_loc
    complex(dp) :: dB_s, dE_s, v_rad
    real(dp)    :: v_sq, T_alpha_erg, A1, A2
    real(dp)    :: mass, charge
    type(profile_values_t) :: pv

    wdot_s = (0.0_dp, 0.0_dp)
    z_cell = x - tetra_physics(ind_tetr)%x1

    s_loc = eval_s_local(ind_tetr, x)
    ! theta_SFL interpolated at the marker position via barycentric weights
    ! over the 4 tetra vertices. Previously we used only the first vertex,
    ! which biased the helical phase factor by up to one tetra in theta and
    ! produced a +(few)e-3 outward shift of the resonance peak in s.
    theta_loc = eval_theta_sfl_local(ind_tetr, x)

    if (coord_system .eq. 1) then
        phi_loc = x(2)
    else
        phi_loc = x(3)
    end if

    B0_loc = tetra_physics(ind_tetr)%bmod1 + sum(tetra_physics(ind_tetr)%gB * z_cell)
    if (B0_loc <= 0.0_dp) return

    call eval_delta_B_s(s_loc, theta_loc, phi_loc, dB_s)
    v_rad = vpar * dB_s
    if (boole_e_perp) then
        call eval_delta_E_s(s_loc, theta_loc, phi_loc, dE_s)
        v_rad = v_rad + clight * dE_s
    end if

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

    ! A_1, A_2 in d/ds form (no ds/dpsi_pol conversion); the contravariant
    ! component delta_B^s in the source closes the chain rule.
    if (charge > 0.0_dp) then
        A1 = pv%dlnn_ds + (charge/T_alpha_erg) * pv%dPhi0_ds - 1.5_dp * pv%dlnTi_ds
        A2 = pv%dlnTi_ds
    else
        A1 = pv%dlnn_ds + (charge/T_alpha_erg) * pv%dPhi0_ds - 1.5_dp * pv%dlnTe_ds
        A2 = pv%dlnTe_ds
    end if

    wdot_s = (-(A1 + A2 * mass * v_sq / (2.0_dp * T_alpha_erg))) &
             * (v_rad / cmplx(B0_loc, 0.0_dp, kind=dp))

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
    use profile_data_mod, only: eval_profiles, profile_values_t, &
                               kim_nu_loaded, eval_kim_nu

    integer,  intent(in) :: n, ind_tetr, species
    real(dp), intent(in) :: x(3), vpar, vperp

    real(dp) :: s_loc, n_loc, T_loc_erg, mass, charge, nu_c, tau_c_loc
    type(profile_values_t) :: pv

    s_loc = eval_s_local(ind_tetr, x)

    if (kim_nu_loaded) then
        ! Use KIM's tabulated nu_e(s) so the regularisation rate is
        ! consistent with the OU step's collisional decorrelation rate.
        nu_c = eval_kim_nu(s_loc)
    else
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
    end if
    if (nu_c <= 0.0_dp) nu_c = 1.0_dp
    tau_c_loc = 1.0_dp / nu_c

    tau_c(n, species)      = tau_c_loc
    nu_r(n, species)       = nu_r_frac * nu_c
    t_reg_on(n, species)   = real(m_collision_times_reg_on, dp) * tau_c_loc

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
! theta_SFL at the marker's cylindrical position x via barycentric
! interpolation over the 4 tetra vertices. Replaces the previous
! "first vertex" lookup which biased the helical phase factor in
! eval_wdot_s by up to one tetra-width in theta and produced a small
! outward shift of the resonance peak in s. Falls back to the first
! vertex value if the linear solve is degenerate.
! ====================================================================
real(dp) function eval_theta_sfl_local(ind_tetr, x) result(theta_loc)

    use tetra_grid_mod, only: tetra_grid, verts_rphiz, verts_sthetaphi
    use tetra_physics_mod, only: coord_system
    use constants, only: pi

    integer,  intent(in) :: ind_tetr
    real(dp), intent(in) :: x(3)

    integer  :: i, vidx, ipiv(4), ierr
    real(dp) :: M(4,4), b(4,1), V_th(4), span

    ! In SFL coordinates (coord_system==2), theta is the native x(2) coordinate.
    if (coord_system == 2) then
        theta_loc = modulo(x(2), 2.0_dp*pi)
        return
    end if

    do i = 1, 4
        vidx = tetra_grid(ind_tetr)%ind_knot(i)
        M(1:3, i) = verts_rphiz(:, vidx)
        M(4,   i) = 1.0_dp
        V_th(i)   = verts_sthetaphi(2, vidx)
    end do
    b(1:3, 1) = x
    b(4,   1) = 1.0_dp

    call dgesv(4, 1, M, 4, ipiv, b, 4, ierr)
    if (ierr /= 0) then
        theta_loc = V_th(1)
        return
    end if

    ! Unwrap theta across the 2pi seam if the tetra straddles it.
    span = maxval(V_th) - minval(V_th)
    if (span > pi) then
        do i = 1, 4
            if (V_th(i) < pi) V_th(i) = V_th(i) + 2.0_dp * pi
        end do
    end if

    theta_loc = sum(b(:, 1) * V_th)
    theta_loc = modulo(theta_loc, 2.0_dp * pi)

end function eval_theta_sfl_local

! ====================================================================
! Total physical volume of the spawn region (the part of the mesh whose
! prisms have their reference s in [s_inner_sample, s_outer_sample]).
! This is the V_W factor that the delta-f current density is missing
! before applying it to the moments (see
! docs/2026-05-06-spawn-volume-normalisation.md).
!
! For each prism we use verts_sthetaphi at the first knot of the prism's
! first tetra as the prism's reference s -- same convention as the
! diagnostic scripts. Sums over all phi slices since markers are placed
! uniformly in phi too.
! ====================================================================
real(dp) function compute_spawn_volume() result(V_W)

    use gorilla_applets_types_mod, only: output
    use tetra_grid_mod, only: tetra_grid, verts_sthetaphi, ntetr

    integer  :: n_prisms, p, knot1
    real(dp) :: s_p

    n_prisms = ntetr / 3
    V_W = 0.0_dp
    do p = 1, n_prisms
        knot1 = tetra_grid(3*p - 2)%ind_knot(1)
        s_p   = verts_sthetaphi(1, knot1)
        if (s_p >= s_inner_sample .and. s_p <= s_outer_sample) then
            V_W = V_W + output%prism_volumes(p)
        end if
    end do

end function compute_spawn_volume

! ====================================================================
! Population filter: classify each marker as trapped or passing using
! conservation of energy and magnetic moment, then mark non-matching
! ones as lost. Trapping condition:
!     pitch^2 < 1 - B_start / B_max(s_marker),
! with B_max(s) found by scanning theta on the marker's flux surface
! via magdata_in_symfluxcoord_ext.
!
! Note on normalisation: filtered-out markers stay in the n_particles
! divisor, so the plotted magnitude is the FILTERED population's
! contribution to the (unfiltered) j_par. Passing-only and trapped-only
! plotted profiles add up to the original.
! ====================================================================
subroutine filter_markers_by_trapping(species)

    use gorilla_applets_types_mod, only: in, start
    use tetra_physics_mod, only: tetra_physics
    use find_tetra_mod, only: find_tetra
    use magdata_in_symfluxcoordinates_mod, only: magdata_in_symfluxcoord_ext
    use constants, only: pi

    integer, intent(in) :: species

    integer, parameter :: n_theta_scan = 100
    integer  :: n, ind_tetr, iface, k, n_kept, n_dropped
    real(dp) :: x(3), z_cell(3), s_loc, B_start, B_max
    real(dp) :: pitch, eps_eff, theta_scan
    real(dp) :: psi_d, q_d, dq_d, sqrtg_d, dbmod_dt_d
    real(dp) :: bmod_at_scan, R_d, dR_ds_d, dR_dt_d, Z_d, dZ_ds_d, dZ_dt_d
    logical  :: is_trapped, keep

    if (trim(trapping_filter_mode) /= 'passing' .and. &
        trim(trapping_filter_mode) /= 'trapped') return

    n_kept    = 0
    n_dropped = 0

    do n = 1, in%num_particles
        if (start%lost(n, species)) cycle

        x = start%x(:, n, species)
        call find_tetra(x, 0.0_dp, 0.0_dp, ind_tetr, iface)
        if (ind_tetr == -1) then
            start%lost(n, species) = .true.
            n_dropped = n_dropped + 1
            cycle
        end if

        s_loc   = eval_s_local(ind_tetr, x)
        z_cell  = x - tetra_physics(ind_tetr)%x1
        B_start = tetra_physics(ind_tetr)%bmod1 + sum(tetra_physics(ind_tetr)%gB * z_cell)

        ! B_max on the flux surface s_loc.
        B_max = 0.0_dp
        do k = 1, n_theta_scan
            theta_scan = 2.0_dp * pi * (k - 1) / n_theta_scan
            psi_d = 0.0_dp
            call magdata_in_symfluxcoord_ext(1, s_loc, psi_d, theta_scan, &
                                             q_d, dq_d, sqrtg_d, bmod_at_scan, dbmod_dt_d, &
                                             R_d, dR_ds_d, dR_dt_d, &
                                             Z_d, dZ_ds_d, dZ_dt_d)
            if (bmod_at_scan > B_max) B_max = bmod_at_scan
        end do

        pitch   = start%pitch(n, species)
        eps_eff = 1.0_dp - B_start / B_max
        is_trapped = (eps_eff > 0.0_dp .and. pitch * pitch < eps_eff)

        if (trim(trapping_filter_mode) == 'passing') then
            keep = .not. is_trapped
        else
            keep = is_trapped
        end if

        if (.not. keep) then
            start%lost(n, species) = .true.
            n_dropped = n_dropped + 1
        else
            n_kept = n_kept + 1
        end if
    end do

    print *, ''
    print '(a, a)', ' [trapping filter] mode = ', trim(trapping_filter_mode)
    print '(a, i0, a, i0)', '   kept: ', n_kept, ', dropped: ', n_dropped

end subroutine filter_markers_by_trapping

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
    use tetra_grid_settings_mod, only: grid_kind
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
                if (coord_system .eq. 1 .and. grid_kind .eq. 5) then
                    ! Analytic-circ: draw directly inside the annulus.
                    call draw_annulus_rphiz_analytic(s_inner_sample, s_outer_sample, u, x)
                    call find_tetra(x, 0.0_dp, 0.0_dp, ind_tetr, iface)
                    if (ind_tetr == -1) cycle
                    start%x(:, n, species) = x
                    n_replaced = n_replaced + 1
                    exit
                else if (coord_system .eq. 1) then
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

! ====================================================================
! Equidistant-in-s spawner. Reconstructs the radial s-edges of the grid
! (matching tetra_grid_settings_mod's logarithmic-extra-rings + linear
! main-grid layout) and places ONE marker at each layer whose midpoint
! falls inside [s_inner_sample, s_outer_sample]. Theta and phi are
! sampled uniformly from the full mesh box, and a position is accepted
! only if its local s lands inside that layer (rejection per layer).
!
! n_particles in the input must be >= the number of layers in the
! window; layers in excess of the slot count are skipped and a warning
! is printed. Unused slots get marked as lost.
! ====================================================================
subroutine spawn_equidistant_in_s(species, n_spawned)

    use gorilla_applets_types_mod, only: in, start, g
    use tetra_grid_settings_mod, only: sfc_s_min, sfc_s_max, grid_size, n_extra_rings, &
                                       grid_kind, R0_analytic_circ, a_analytic_circ
    use tetra_physics_mod, only: coord_system
    use find_tetra_mod, only: find_tetra
    use magdata_in_symfluxcoordinates_mod, only: magdata_in_symfluxcoord_ext
    use constants, only: pi

    integer, intent(in)  :: species
    integer, intent(out) :: n_spawned

    ! Direct (s, theta, phi) -> (R, phi, Z) sampler via the equilibrium's
    ! symflux mapping. No rejection, so a small retry budget covers the
    ! occasional find_tetra miss at the very edge of the mesh.
    integer, parameter   :: max_tries = 20
    integer  :: nr, n1, k_layer, ind_tetr, iface, n_tries, k
    integer  :: n_layers_window, k_per_layer, k_in_layer
    integer  :: j_swap, tmp_int, phi_bin_idx
    integer, allocatable :: phi_bin_perm(:)
    real(dp) :: s_mid, s_low_lyr, s_high_lyr, s_target, x(3), u(3), u_scalar
    real(dp) :: theta_loc, phi_loc, R_loc, Z_loc
    real(dp) :: s_second_ring
    real(dp), allocatable :: s_edges(:)
    logical  :: placed
    ! Analytic circular tokamak (grid_kind=5) closed-form (s,theta)->(R,Z) map.
    real(dp) :: s_edge_ac, rho_loc
    ! Dummy outputs from magdata_in_symfluxcoord_ext (we only need R, Z).
    real(dp) :: psi_d, q_d, dq_ds_d, sqrtg_d, bmod_d, dbmod_dt_d
    real(dp) :: dR_ds_d, dR_dt_d, dZ_ds_d, dZ_dt_d

    nr = grid_size(1)
    n1 = nr - n_extra_rings

    allocate(s_edges(nr + 1))
    s_edges(1) = sfc_s_min
    if (n_extra_rings > 0) then
        s_second_ring = sfc_s_min + (sfc_s_max - sfc_s_min) / dble(n1)
        do k = 1, n_extra_rings
            s_edges(k + 1) = exp(log(sfc_s_min) &
                + dble(k) * (log(s_second_ring) - log(sfc_s_min)) / dble(n_extra_rings + 1))
        end do
    end if
    do k = 1, n1
        s_edges(n_extra_rings + 1 + k) = sfc_s_min + dble(k) * (sfc_s_max - sfc_s_min) / dble(n1)
    end do

    ! Pre-pass: count layers whose midpoint falls in the spawn window.
    ! Sets the multiplicity per layer when n_particles > n_layers.
    n_layers_window = 0
    do k_layer = 1, nr
        s_mid = 0.5_dp * (s_edges(k_layer) + s_edges(k_layer + 1))
        if (s_mid >= s_inner_sample .and. s_mid <= s_outer_sample) then
            n_layers_window = n_layers_window + 1
        end if
    end do

    if (n_layers_window == 0) then
        print *, 'ERROR (spawn_equidistant_in_s): no layer midpoints fall in spawn window.'
        deallocate(s_edges)
        n_spawned = 0
        return
    end if

    ! Markers per layer: ceil(n_particles / n_layers_window). The very last
    ! layer may get fewer than this if n_particles isn't an exact multiple
    ! (we cap at in%num_particles).
    k_per_layer = (in%num_particles + n_layers_window - 1) / n_layers_window

    print *, ''
    print *, 'Equidistant-in-s spawn enabled.'
    print '(a, i0, a, i0)', '   Total radial layers: ', nr, ', main-grid layers (n1): ', n1
    print '(a, f8.5, a, f8.5, a)', '   Spawn window: [', s_inner_sample, &
                                    ', ', s_outer_sample, ']'
    print '(a, i0, a, i0, a, i0)', '   Layers in window: ', n_layers_window, &
                                    ', markers per layer: ', k_per_layer, &
                                    ' (total slots: ', in%num_particles, ')'

    n_spawned = 0

    if (boole_stratify_phi) allocate(phi_bin_perm(k_per_layer))

    do k_layer = 1, nr
        s_low_lyr  = s_edges(k_layer)
        s_high_lyr = s_edges(k_layer + 1)
        s_mid      = 0.5_dp * (s_low_lyr + s_high_lyr)

        if (s_mid < s_inner_sample .or. s_mid > s_outer_sample) cycle

        if (boole_stratify_phi) then
            ! Fresh Fisher-Yates shuffle of bin indices for this layer.
            do k = 1, k_per_layer
                phi_bin_perm(k) = k
            end do
            do k = k_per_layer, 2, -1
                call random_number(u_scalar)
                j_swap = 1 + int(u_scalar * dble(k))
                if (j_swap > k) j_swap = k
                tmp_int           = phi_bin_perm(k)
                phi_bin_perm(k)   = phi_bin_perm(j_swap)
                phi_bin_perm(j_swap) = tmp_int
            end do
        end if

        do k_in_layer = 1, k_per_layer
            if (n_spawned >= in%num_particles) exit

            placed = .false.
            do n_tries = 1, max_tries
                call random_number(u)
                ! Direct draw inside this layer: s uniform within layer
                ! edges, theta and phi uniform in [0, 2 pi]. With
                ! boole_stratify_theta on, marker k of K in this layer
                ! draws theta from bin [(k-1)/K, k/K)*2pi instead.
                s_target  = s_low_lyr + u(1) * (s_high_lyr - s_low_lyr)
                if (boole_stratify_theta) then
                    theta_loc = (dble(k_in_layer - 1) + u(2)) &
                              * 2.0_dp * pi / dble(k_per_layer)
                else
                    theta_loc = u(2) * 2.0_dp * pi
                end if
                if (boole_stratify_phi) then
                    phi_bin_idx = phi_bin_perm(k_in_layer)
                    phi_loc = (dble(phi_bin_idx - 1) + u(3)) &
                            * 2.0_dp * pi / dble(k_per_layer)
                else
                    phi_loc = u(3) * 2.0_dp * pi
                end if

                if (coord_system == 1) then
                    if (grid_kind == 5) then
                        ! Analytic circular tokamak: closed-form (s,theta_SFL)->(R,Z).
                        ! magdata_in_symfluxcoord_ext is EFIT-only (its splines are
                        ! never loaded for grid_kind=5), so we map directly with the
                        ! same geometry the grid builder uses (tetra_grid_mod case(5)):
                        !   s_edge = R0 - sqrt(R0^2 - a^2);  rho = sqrt(R0^2 - (R0 - s*s_edge)^2)
                        s_edge_ac = R0_analytic_circ &
                                  - sqrt(R0_analytic_circ**2 - a_analytic_circ**2)
                        rho_loc = sqrt(R0_analytic_circ**2 &
                                     - (R0_analytic_circ - s_target*s_edge_ac)**2)
                        R_loc = R0_analytic_circ + rho_loc * cos(theta_loc)
                        Z_loc = rho_loc * sin(theta_loc)
                    else
                        ! Symflux -> cylindrical via the equilibrium mapping.
                        ! psi_d is overwritten as output when inp_label=1.
                        psi_d = 0.0_dp
                        call magdata_in_symfluxcoord_ext(1, s_target, psi_d, theta_loc, &
                                                         q_d, dq_ds_d, sqrtg_d, bmod_d, dbmod_dt_d, &
                                                         R_loc, dR_ds_d, dR_dt_d, &
                                                         Z_loc, dZ_ds_d, dZ_dt_d)
                    end if
                    x(1) = R_loc
                    x(2) = phi_loc
                    x(3) = Z_loc
                else
                    x(1) = s_target
                    x(2) = theta_loc
                    x(3) = phi_loc
                end if

                call find_tetra(x, 0.0_dp, 0.0_dp, ind_tetr, iface)
                if (ind_tetr == -1) cycle  ! at the mesh edge; retry with a fresh draw

                n_spawned = n_spawned + 1
                start%x(:, n_spawned, species) = x
                start%lost(n_spawned, species) = .false.
                placed = .true.
                exit
            end do
            if (.not. placed) then
                print '(a, i0, a, i0, a, f7.4, a)', &
                    '   WARNING: could not place marker ', k_in_layer, &
                    ' for layer ', k_layer, ' (s_mid = ', s_mid, ')'
            end if
        end do
    end do

    do k = n_spawned + 1, in%num_particles
        start%lost(k, species) = .true.
    end do

    deallocate(s_edges)
    if (allocated(phi_bin_perm)) deallocate(phi_bin_perm)

    print '(a, i0, a, i0, a)', '   Placed ', n_spawned, &
                                ' markers (', in%num_particles, ' slots available)'
    print *, ''

end subroutine spawn_equidistant_in_s

! ====================================================================
! Analytic circular tokamak (grid_kind=5) direct draw of a cylindrical
! (R, phi, Z) position inside the resonance annulus s in [s_lo, s_hi].
! Given three uniform deviates u(1:3), maps (s,theta,phi) -> (R,phi,Z) via
! the closed-form geometry. rho is drawn with pdf ~ rho (rho^2 uniform) so
! the sample is uniform in poloidal AREA within the ring -- identical to the
! distribution the (R,Z)-box rejection sampler produced, but with ~100%
! acceptance instead of ~1% (the ring is a thin sliver of the full poloidal
! bounding box). This restricts respawned markers to the thin resonance
! annulus and removes the marker losses caused by rejection-sampling misses.
! ====================================================================
subroutine draw_annulus_rphiz_analytic(s_lo, s_hi, u, x)

    use tetra_grid_settings_mod, only: R0_analytic_circ, a_analytic_circ, n_field_periods
    use constants, only: pi

    real(dp), intent(in)  :: s_lo, s_hi, u(3)
    real(dp), intent(out) :: x(3)
    real(dp) :: s_edge, rho_lo, rho_hi, rho_s, theta_s

    s_edge = R0_analytic_circ - sqrt(R0_analytic_circ**2 - a_analytic_circ**2)
    rho_lo = sqrt(R0_analytic_circ**2 - (R0_analytic_circ - s_lo*s_edge)**2)
    rho_hi = sqrt(R0_analytic_circ**2 - (R0_analytic_circ - s_hi*s_edge)**2)
    rho_s   = sqrt(rho_lo**2 + u(1)*(rho_hi**2 - rho_lo**2))   ! uniform in area
    theta_s = u(2) * 2.0_dp * pi                     ! poloidal: always full [0,2*pi]
    x(1) = R0_analytic_circ + rho_s * cos(theta_s)   ! R
    ! Toroidal: sample only the modelled field period [0, 2*pi/n_field_periods]
    ! so respawned/spawned markers land inside the wedge mesh (n_field_periods=N).
    ! n_field_periods = 1 -> full torus -> exact no-op.
    x(2) = u(3) * 2.0_dp * pi / dble(n_field_periods) ! phi
    x(3) = rho_s * sin(theta_s)                      ! Z

end subroutine draw_annulus_rphiz_analytic

! ====================================================================
! Single-particle respawn within the bias window intersected with the
! computational-domain cut. Used to redraw start%x(:,n,species) when a
! marker exits the domain mid-trace, so it can be re-traced from a fresh
! valid position. Thread-safe (only writes to its own particle's start
! position; find_tetra reads the immutable mesh; random_number is per
! thread under gfortran/OpenMP).
! ====================================================================
subroutine respawn_particle_in_window(n, species, success)

    use gorilla_applets_types_mod, only: start, g
    use tetra_physics_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind
    use find_tetra_mod, only: find_tetra
    use constants, only: pi

    integer, intent(in)  :: n, species
    logical, intent(out) :: success

    integer, parameter :: max_tries = 1000
    integer  :: ind_tetr, iface, n_tries
    real(dp) :: x(3), s_loc, u(3), x_amin, x_amax, x_cmin, x_cmax
    real(dp) :: s_lo, s_hi
    real(dp), parameter :: eps_cut = 1.0e-3_dp

    s_lo = s_inner_sample
    s_hi = s_outer_sample
    if (s_outer_cut < 1.0_dp) s_hi = min(s_hi, s_outer_cut - eps_cut)

    success = .false.
    if (s_hi <= s_lo) return

    x_amin = g%amin; x_amax = g%amax
    x_cmin = g%cmin; x_cmax = g%cmax

    do n_tries = 1, max_tries
        call random_number(u)
        if (coord_system .eq. 1 .and. grid_kind .eq. 5) then
            ! Analytic-circ: draw directly inside the annulus (s in [s_lo,s_hi]
            ! guaranteed by construction), so no s-rejection is needed.
            call draw_annulus_rphiz_analytic(s_lo, s_hi, u, x)
            call find_tetra(x, 0.0_dp, 0.0_dp, ind_tetr, iface)
            if (ind_tetr == -1) cycle   ! rare find_tetra miss at a cell edge
            start%x(:, n, species) = x
            success = .true.
            return
        else if (coord_system .eq. 1) then
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
        if (s_loc >= s_lo .and. s_loc <= s_hi) then
            start%x(:, n, species) = x
            success = .true.
            return
        end if
    end do

end subroutine respawn_particle_in_window

! ====================================================================
! Diagnostic: dump per-particle starting positions to a file. Writes
! columns
!   n  x1  x2  x3  s  lost
! where (x1, x2, x3) are in the active coord_system (cylindrical:
! (R, phi, Z); flux: (s, theta, phi)). s is evaluated via find_tetra
! at the same position; -1 marks particles outside the mesh.
! ====================================================================
subroutine dump_start_positions(filename, species)

    use gorilla_applets_types_mod, only: in, start
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: tetra_grid, verts_sthetaphi
    use find_tetra_mod, only: find_tetra

    character(len=*), intent(in) :: filename
    integer, intent(in)          :: species

    integer  :: u, n, ind_tetr, iface, lost_int
    real(dp) :: x(3), s_loc, theta_loc, z_cell(3), B_start, pitch

    open(newunit=u, file=trim(filename), status='unknown', action='write')
    write(u, '(a)') '# n  x1  x2  x3  s  theta_SFL  B_start  pitch  lost'
    do n = 1, in%num_particles
        x = start%x(:, n, species)
        s_loc      = -1.0_dp
        theta_loc  = -1.0_dp
        B_start    = -1.0_dp
        pitch      =  start%pitch(n, species)
        if (.not. start%lost(n, species)) then
            call find_tetra(x, 0.0_dp, 0.0_dp, ind_tetr, iface)
            if (ind_tetr /= -1) then
                s_loc = eval_s_local(ind_tetr, x)
                ! theta_SFL via barycentric interpolation at the marker's
                ! cylindrical position over the 4 tetra vertices (same
                ! convention as eval_wdot_s).
                theta_loc = eval_theta_sfl_local(ind_tetr, x)
                z_cell = x - tetra_physics(ind_tetr)%x1
                B_start = tetra_physics(ind_tetr)%bmod1 + sum(tetra_physics(ind_tetr)%gB * z_cell)
            end if
        end if
        if (start%lost(n, species)) then
            lost_int = 1
        else
            lost_int = 0
        end if
        write(u, '(i6, 7es16.8, i4)') n, x(1), x(2), x(3), s_loc, theta_loc, B_start, pitch, lost_int
    end do
    close(u)

end subroutine dump_start_positions

! ====================================================================
! Velocity-Jacobian-weighted Maxwellian energy PDF for the 1D sampler.
! p(E) dE = (4 pi v^2 / m) * f_M(v) dE  ∝  sqrt(E) * exp(-E/T)
! T is taken from module-private T_sample_eV, which the caller (the
! 'maxwellian' branch of set_rest_of_start_type_rmp) sets from
! in%energy_eV before calling init_distribution_1d.  Matches
! pdf_interface(x(:)) -> val so it can be passed to init_distribution_1d.
! ====================================================================
function pdf_maxwellian_energy(x) result(val)
    real(dp), intent(in) :: x(:)
    real(dp) :: val, E

    E = x(1)
    if (E <= 0.0_dp .or. T_sample_eV <= 0.0_dp) then
        val = 0.0_dp
        return
    endif
    val = sqrt(E) * exp(-E / T_sample_eV)
end function pdf_maxwellian_energy

end module utils_rmp_response_currents_mod
