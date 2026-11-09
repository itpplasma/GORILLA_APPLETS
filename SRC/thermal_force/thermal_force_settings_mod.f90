module thermal_force_settings_mod

! Input of the thermal-force applet (i_option = 16), file thermal_force.inp,
! namelist thermal_force_nml. Units: CGS with temperatures in eV, masses in
! amu, as in the names. See docs/plans/2026-11-05-homma-thermal-force-design.md.

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    integer, parameter, public :: max_output_times = 256

    type, public :: thermal_force_input_t
        ! 'velocity': homogeneous plasma, Cartesian velocities with Lorentz rotation
        ! 'slab': uniform B, linear temperature profile, ions move in space
        character(len=32) :: driver = 'velocity'
        ! slab driver: 'guiding_centre' (GORILLA model, homma_gc_step) or
        ! 'full_orbit' (reference only, not GORILLA physics: Homma's full-orbit
        ! model, exact helix between Cartesian kicks at the particle; needs
        ! Omega_a dt <~ 1, otherwise the kicks alias with the gyration)
        character(len=32) :: orbit_model = 'guiding_centre'
        ! test species
        real(dp) :: test_mass_amu = 183.84_dp
        real(dp) :: test_charge_number = 3.0_dp
        ! background ions (one species, no electrons)
        real(dp) :: background_mass_amu = 1.007276466621_dp
        real(dp) :: background_charge_number = 1.0_dp
        real(dp) :: density_cm3 = 1.0e14_dp
        real(dp) :: temperature_ev = 50.0_dp
        real(dp) :: flow_cm_s(3) = 0.0_dp
        real(dp) :: grad_temperature_ev_cm(3) = 0.0_dp
        real(dp) :: coulomb_log = 15.0_dp
        character(len=32) :: heat_flux_model = 'homma2014'          ! or 'none'
        character(len=32) :: collision_background_at = 'gyro_orbit' ! or 'guiding_centre'
        ! magnetic field (uniform in the velocity and slab drivers)
        real(dp) :: b_gauss(3) = [0.0_dp, 0.0_dp, 1.0e4_dp]
        ! initial test distribution: 'maxwellian' or 'drifting_bimaxwellian'
        ! (parallel and perpendicular to b), drift vector in the lab frame
        character(len=32) :: test_distribution = 'maxwellian'
        real(dp) :: test_temperature_ev = 50.0_dp
        real(dp) :: test_temperature_par_ev = 50.0_dp
        real(dp) :: test_temperature_perp_ev = 50.0_dp
        real(dp) :: test_drift_cm_s(3) = 0.0_dp
        real(dp) :: start_x_cm(3) = 0.0_dp
        ! collision step in units of tau_s at v = sqrt(3 T_b/m_a)
        real(dp) :: dt_coll_over_tau_s = 0.01_dp
        ! output times in 'tau_s' or 'gyro' (Omega_a t) units, increasing;
        ! unused entries stay negative
        character(len=32) :: output_time_unit = 'tau_s'
        real(dp) :: output_times(max_output_times) = -1.0_dp
        ! Gyration of the velocity between kicks (velocity driver)
        logical :: boole_gyration = .true.
        ! Also run grad T = 0 with the same random numbers; output the paired difference
        logical :: boole_control_run = .true.
        integer :: n_particles = 10000
        integer :: seed = 1
        character(len=256) :: output_file = 'thermal_force_velocity.csv'
        ! slab driver: positions of the first n_trajectories ions (gradient run)
        ! every trajectory_dt_gyro gyro periods (full orbit) or collision step
        integer :: n_trajectories = 0
        real(dp) :: trajectory_dt_gyro = 0.25_dp
        character(len=256) :: trajectory_file = 'thermal_force_trajectories.csv'
    end type thermal_force_input_t

    public :: load_thermal_force_inp

contains

subroutine load_thermal_force_inp(inp)

    type(thermal_force_input_t), intent(out) :: inp
    character(len=32) :: driver, orbit_model, heat_flux_model, collision_background_at, test_distribution, output_time_unit
    real(dp) :: test_mass_amu, test_charge_number, background_mass_amu, background_charge_number
    real(dp) :: density_cm3, temperature_ev, flow_cm_s(3), grad_temperature_ev_cm(3), coulomb_log, b_gauss(3)
    real(dp) :: test_temperature_ev, test_temperature_par_ev, test_temperature_perp_ev, test_drift_cm_s(3)
    real(dp) :: start_x_cm(3), dt_coll_over_tau_s, output_times(max_output_times)
    logical :: boole_gyration, boole_control_run
    integer :: n_particles, seed, unit, n_trajectories
    character(len=256) :: output_file, trajectory_file
    real(dp) :: trajectory_dt_gyro

    namelist /thermal_force_nml/ driver, test_mass_amu, test_charge_number, background_mass_amu, &
        background_charge_number, density_cm3, temperature_ev, flow_cm_s, grad_temperature_ev_cm, &
        coulomb_log, heat_flux_model, collision_background_at, b_gauss, test_distribution, &
        test_temperature_ev, test_temperature_par_ev, test_temperature_perp_ev, test_drift_cm_s, &
        start_x_cm, dt_coll_over_tau_s, output_time_unit, output_times, boole_gyration, &
        boole_control_run, n_particles, seed, output_file, orbit_model, n_trajectories, &
        trajectory_dt_gyro, trajectory_file

    ! defaults from the type
    driver = inp%driver; test_mass_amu = inp%test_mass_amu; test_charge_number = inp%test_charge_number
    background_mass_amu = inp%background_mass_amu; background_charge_number = inp%background_charge_number
    density_cm3 = inp%density_cm3; temperature_ev = inp%temperature_ev; flow_cm_s = inp%flow_cm_s
    grad_temperature_ev_cm = inp%grad_temperature_ev_cm; coulomb_log = inp%coulomb_log
    heat_flux_model = inp%heat_flux_model; collision_background_at = inp%collision_background_at
    b_gauss = inp%b_gauss; test_distribution = inp%test_distribution
    test_temperature_ev = inp%test_temperature_ev; test_temperature_par_ev = inp%test_temperature_par_ev
    test_temperature_perp_ev = inp%test_temperature_perp_ev; test_drift_cm_s = inp%test_drift_cm_s
    start_x_cm = inp%start_x_cm; dt_coll_over_tau_s = inp%dt_coll_over_tau_s
    output_time_unit = inp%output_time_unit; output_times = inp%output_times
    boole_gyration = inp%boole_gyration; boole_control_run = inp%boole_control_run
    n_particles = inp%n_particles; seed = inp%seed; output_file = inp%output_file
    orbit_model = inp%orbit_model; n_trajectories = inp%n_trajectories
    trajectory_dt_gyro = inp%trajectory_dt_gyro; trajectory_file = inp%trajectory_file

    open(newunit=unit, file='thermal_force.inp', status='old', action='read')
    read(unit, nml=thermal_force_nml)
    close(unit)

    inp%driver = driver; inp%test_mass_amu = test_mass_amu; inp%test_charge_number = test_charge_number
    inp%background_mass_amu = background_mass_amu; inp%background_charge_number = background_charge_number
    inp%density_cm3 = density_cm3; inp%temperature_ev = temperature_ev; inp%flow_cm_s = flow_cm_s
    inp%grad_temperature_ev_cm = grad_temperature_ev_cm; inp%coulomb_log = coulomb_log
    inp%heat_flux_model = heat_flux_model; inp%collision_background_at = collision_background_at
    inp%b_gauss = b_gauss; inp%test_distribution = test_distribution
    inp%test_temperature_ev = test_temperature_ev; inp%test_temperature_par_ev = test_temperature_par_ev
    inp%test_temperature_perp_ev = test_temperature_perp_ev; inp%test_drift_cm_s = test_drift_cm_s
    inp%start_x_cm = start_x_cm; inp%dt_coll_over_tau_s = dt_coll_over_tau_s
    inp%output_time_unit = output_time_unit; inp%output_times = output_times
    inp%boole_gyration = boole_gyration; inp%boole_control_run = boole_control_run
    inp%n_particles = n_particles; inp%seed = seed; inp%output_file = output_file
    inp%orbit_model = orbit_model; inp%n_trajectories = n_trajectories
    inp%trajectory_dt_gyro = trajectory_dt_gyro; inp%trajectory_file = trajectory_file

    call check_input(inp)
    print *, 'GORILLA_APPLETS: Loaded input data from thermal_force.inp'

end subroutine load_thermal_force_inp

subroutine check_input(inp)

    type(thermal_force_input_t), intent(in) :: inp
    integer :: n_out

    select case (trim(inp%driver))
    case ('velocity', 'slab')
    case default
        call input_error('driver must be ''velocity'' or ''slab'' (tokamak driver follows)')
    end select
    select case (trim(inp%orbit_model))
    case ('guiding_centre', 'full_orbit')
    case default
        call input_error('orbit_model must be ''guiding_centre'' or ''full_orbit''')
    end select
    if (inp%n_trajectories < 0) call input_error('n_trajectories must not be negative')
    if (inp%trajectory_dt_gyro <= 0.0_dp) call input_error('trajectory_dt_gyro must be positive')
    select case (trim(inp%heat_flux_model))
    case ('homma2014', 'none')
    case default
        call input_error('heat_flux_model must be ''homma2014'' or ''none''')
    end select
    select case (trim(inp%collision_background_at))
    case ('gyro_orbit', 'guiding_centre')
    case default
        call input_error('collision_background_at must be ''gyro_orbit'' or ''guiding_centre''')
    end select
    select case (trim(inp%test_distribution))
    case ('maxwellian', 'drifting_bimaxwellian')
    case default
        call input_error('test_distribution must be ''maxwellian'' or ''drifting_bimaxwellian''')
    end select
    select case (trim(inp%output_time_unit))
    case ('tau_s', 'gyro')
    case default
        call input_error('output_time_unit must be ''tau_s'' or ''gyro''')
    end select
    n_out = count(inp%output_times >= 0.0_dp)
    if (n_out == 0) call input_error('output_times: at least one time >= 0 required')
    if (any(inp%output_times(1:n_out) < 0.0_dp)) call input_error('output_times: negative entry inside the list')
    if (n_out > 1) then
        if (any(inp%output_times(2:n_out) <= inp%output_times(1:n_out-1))) &
            call input_error('output_times must increase')
    end if
    if (inp%dt_coll_over_tau_s <= 0.0_dp) call input_error('dt_coll_over_tau_s must be positive')
    if (inp%n_particles < 2) call input_error('n_particles must be at least 2')
    if (norm2(inp%b_gauss) <= 0.0_dp) call input_error('b_gauss must not vanish')

end subroutine check_input

subroutine input_error(msg)

    character(len=*), intent(in) :: msg

    print *, 'ERROR in thermal_force.inp: ', msg
    error stop 1

end subroutine input_error

end module thermal_force_settings_mod
