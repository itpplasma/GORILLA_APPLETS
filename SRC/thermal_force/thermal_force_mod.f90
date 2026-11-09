module thermal_force_mod

! Thermal-force applet (i_option = 16): Monte Carlo test ions in a background
! with the Homma-Hatayama collision operator (homma_collision_mod).
!
! Driver 'velocity': homogeneous plasma, no spatial motion. Each test ion
! carries its Cartesian velocity, receives Euler-Maruyama collision kicks with
! the background at start_x_cm and gyrates exactly between kicks. This is the
! homogeneous limit of a full orbit and gives the benchmark cases B04
! (relaxation), B05 and B06 (initial thermal-force acceleration).
!
! Driver 'slab': uniform B, temperature linear in space (B07 temperature
! screening, Homma 2014 Figs. 1-4). orbit_model 'guiding_centre' is the
! GORILLA model (homma_gc_step: random gyrophase, kick, shift
! dX = -b x dv/Omega_a, then streaming along b). 'full_orbit' is a reference
! for checking the guiding-centre reduction, not GORILLA physics: Homma's
! model, Cartesian kick with the background at the particle, then the exact
! helix until the next collision (requires Omega_a dt <~ 1).
!
! Driver 'tokamak': GORILLA guiding-centre orbits on its tetrahedral mesh in
! cylindrical coordinates (coord_system = 1), e.g. the analytic circular
! tokamak (grid_kind = 5). Collisions as in the slab guiding-centre model, at
! the field direction and |B| of the mesh; the guiding-centre shift moves the
! ion through the mesh with displace_by_straight_line; between collisions
! GORILLA's pusher integrates the orbit. The background temperature is linear
! in R and Z, T = T_0 + grad T . (x - start_x_cm).
!
! With boole_control_run the same ions are also followed with grad T = 0 and
! the same random numbers; the paired difference isolates the thermal force.
! boole_collision_shift = .false. leaves the guiding centre in place at a
! collision (velocity-only collisions, GORILLA's traditional scope).
!
! Units inside: CGS (GORILLA constants_mod), temperatures as energies.
! Design: docs/plans/2026-11-05-homma-thermal-force-design.md

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use homma_collision_mod, only: homma_units_t, homma_linear_profile_t, homma_local_t, &
                                   homma_local_at, homma_kick, homma_perp_basis, homma_gc_step, &
                                   homma_at_guiding_centre, homma_at_gyro_orbit
    use thermal_force_settings_mod, only: thermal_force_input_t, load_thermal_force_inp
    use thermal_force_rng_mod, only: rng_t, rng_seed, rng_normal3, rng_uniform

    implicit none

    private

    real(dp), parameter :: amu_g = 1.66053906660e-24_dp     ! CODATA 2018
    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

    ! Velocity moments accumulated per output time
    integer, parameter :: q_vx = 1, q_vy = 2, q_vz = 3, q_vpar = 4, q_vpar2 = 5, q_vperp2 = 6, q_v4 = 7
    integer, parameter :: n_q = 7
    integer, parameter :: run_gradient = 1, run_control = 2, run_difference = 3
    ! Slab moments: position relative to start_x_cm, kinetic energy
    integer, parameter :: s_x = 1, s_y = 2, s_z = 3, s_energy = 4
    integer, parameter :: n_s = 4

    type :: setup_t
        type(homma_units_t) :: units
        real(dp) :: m_a, q_a, omega_a, tau_s, dt
        real(dp) :: b_field(3), b_hat(3), e1(3), e2(3)
        type(homma_linear_profile_t) :: prof_gradient, prof_control
        real(dp), allocatable :: t_out(:)
    end type setup_t

    public :: calc_thermal_force

contains

subroutine calc_thermal_force()

    type(thermal_force_input_t) :: inp
    type(setup_t) :: su

    call load_thermal_force_inp(inp)
    if (trim(inp%driver) == 'tokamak') call initialize_gorilla_for_test_ion(inp)
    call build_setup(inp, su)

    print '(a,es12.5,a,es12.5,a,es12.5,a)', ' thermal force: tau_s = ', su%tau_s, ' s, Omega_a = ', su%omega_a, &
        ' 1/s, collision step = ', su%dt, ' s'

    select case (trim(inp%driver))
    case ('velocity')
        call run_velocity_driver(inp, su)
    case ('slab')
        call run_slab_driver(inp, su)
    case ('tokamak')
        call run_tokamak_driver(inp, su)
    end select

end subroutine calc_thermal_force

subroutine build_setup(inp, su)

    use constants, only: clight, echarge, ev2erg

    type(thermal_force_input_t), intent(in) :: inp
    type(setup_t), intent(out) :: su
    real(dp) :: t_b, m_b, q_b, v, x, mu
    integer :: n_out

    su%units%k_coulomb = 1.0_dp
    su%units%k_lorentz = 1.0_dp/clight

    su%m_a = inp%test_mass_amu*amu_g
    su%q_a = inp%test_charge_number*echarge
    m_b = inp%background_mass_amu*amu_g
    q_b = inp%background_charge_number*echarge
    t_b = inp%temperature_ev*ev2erg

    su%prof_gradient%bg0%mass = m_b
    su%prof_gradient%bg0%charge = q_b
    su%prof_gradient%bg0%density = inp%density_cm3
    su%prof_gradient%bg0%temperature = t_b
    su%prof_gradient%bg0%flow = inp%flow_cm_s
    su%prof_gradient%bg0%coulomb_log = inp%coulomb_log
    su%prof_gradient%x0 = inp%start_x_cm
    su%prof_gradient%grad_t = inp%grad_temperature_ev_cm*ev2erg
    su%prof_gradient%heat_flux = trim(inp%heat_flux_model) == 'homma2014'
    su%prof_control = su%prof_gradient
    su%prof_control%grad_t = 0.0_dp

    su%b_field = inp%b_gauss
    if (trim(inp%driver) == 'tokamak') then
        call tokamak_start_field(inp, su%b_field)
        ! the profile lives in the local physical frame (R, 0, Z), see tokamak_ion
        su%prof_gradient%x0 = [inp%start_x_cm(1), 0.0_dp, inp%start_x_cm(3)]
    end if
    su%b_hat = su%b_field/norm2(su%b_field)
    call homma_perp_basis(su%b_hat, su%e1, su%e2)
    su%omega_a = su%units%k_lorentz*su%q_a*norm2(su%b_field)/su%m_a

    ! Slowing-down time (Homma 2014 Sec. 3.2) at v = sqrt(3 T_b/m_a):
    ! tau_s = m_a^2 v^3/(4 pi k_C^2 (1 + m_a/m_b) mu(x) q_a^2 q_b^2 n lnL),
    ! x = m_b v^2/(2 T_b), mu(x) = erf(sqrt(x)) - 2 sqrt(x/pi) exp(-x)
    v = sqrt(3.0_dp*t_b/su%m_a)
    x = m_b*v**2/(2.0_dp*t_b)
    mu = erf(sqrt(x)) - 2.0_dp*sqrt(x/pi)*exp(-x)
    su%tau_s = su%m_a**2*v**3/(4.0_dp*pi*su%units%k_coulomb**2*(1.0_dp + su%m_a/m_b)*mu &
             *su%q_a**2*q_b**2*inp%density_cm3*inp%coulomb_log)
    su%dt = inp%dt_coll_over_tau_s*su%tau_s

    n_out = count(inp%output_times >= 0.0_dp)
    select case (trim(inp%output_time_unit))
    case ('tau_s')
        su%t_out = inp%output_times(1:n_out)*su%tau_s
    case ('gyro')
        su%t_out = inp%output_times(1:n_out)/abs(su%omega_a)
    end select

end subroutine build_setup

subroutine run_velocity_driver(inp, su)

    use constants, only: ev2erg

    type(thermal_force_input_t), intent(in) :: inp
    type(setup_t), intent(in) :: su
    type(homma_local_t) :: loc_gradient(1), loc_control(1)
    type(rng_t) :: rng
    real(dp), allocatable :: sums(:,:,:,:)
    real(dp) :: v_grad(3), v_ctrl(3), xi(3), dv(3), t, h, q_grad(n_q), q_ctrl(n_q)
    integer(int64) :: fails(2)
    integer :: ip, k, n_out, ierr

    n_out = size(su%t_out)
    allocate(sums(2, n_q, 3, n_out))
    sums = 0.0_dp
    fails = 0

    call homma_local_at(su%units, su%m_a, su%q_a, su%prof_gradient, inp%start_x_cm, su%b_field, &
                        loc_gradient(1), ierr)
    if (ierr /= 0) error stop 'thermal force: background temperature not positive at start_x_cm'
    call homma_local_at(su%units, su%m_a, su%q_a, su%prof_control, inp%start_x_cm, su%b_field, &
                        loc_control(1), ierr)

    !$omp parallel do schedule(dynamic, 256) default(shared) &
    !$omp private(ip, rng, v_grad, v_ctrl, xi, dv, t, h, k, ierr, q_grad, q_ctrl) &
    !$omp reduction(+:sums, fails)
    do ip = 1, inp%n_particles
        call rng_seed(rng, int(inp%seed, int64), int(ip, int64))
        call sample_initial_velocity(inp, su, rng, v_grad)
        v_ctrl = v_grad
        t = 0.0_dp
        do k = 1, n_out
            do while (t < su%t_out(k))
                h = min(su%dt, su%t_out(k) - t)
                call rng_normal3(rng, xi)
                call homma_kick(loc_gradient, v_grad, h, xi, dv, ierr)
                if (ierr /= 0) fails(1) = fails(1) + 1
                v_grad = v_grad + dv
                if (inp%boole_control_run) then
                    call homma_kick(loc_control, v_ctrl, h, xi, dv, ierr)
                    if (ierr /= 0) fails(2) = fails(2) + 1
                    v_ctrl = v_ctrl + dv
                end if
                if (inp%boole_gyration) then
                    call gyrate(su, -su%omega_a*h, v_grad)
                    call gyrate(su, -su%omega_a*h, v_ctrl)
                end if
                if (h < su%dt) then
                    t = su%t_out(k)
                else
                    t = t + h
                end if
            end do
            q_grad = moments(su, v_grad)
            q_ctrl = moments(su, v_ctrl)
            sums(1,:,run_gradient,k) = sums(1,:,run_gradient,k) + q_grad
            sums(2,:,run_gradient,k) = sums(2,:,run_gradient,k) + q_grad**2
            sums(1,:,run_control,k) = sums(1,:,run_control,k) + q_ctrl
            sums(2,:,run_control,k) = sums(2,:,run_control,k) + q_ctrl**2
            sums(1,:,run_difference,k) = sums(1,:,run_difference,k) + (q_grad - q_ctrl)
            sums(2,:,run_difference,k) = sums(2,:,run_difference,k) + (q_grad - q_ctrl)**2
        end do
    end do
    !$omp end parallel do

    if (any(fails > 0)) print '(a,2i12)', ' thermal force: kernel failures (gradient, control) = ', fails
    call write_velocity_output(inp, su, sums, fails, ev2erg)

end subroutine run_velocity_driver

subroutine run_slab_driver(inp, su)

    ! Collisions on the fixed grid t_n = n dt (kick at t_n), exact free motion
    ! in between. Output times sample the free motion and do not shorten the
    ! collision steps.

    type(thermal_force_input_t), intent(in) :: inp
    type(setup_t), intent(in) :: su
    real(dp), allocatable :: sums(:,:,:,:)
    integer(int64) :: fails(2)
    integer :: ip, traj_unit

    allocate(sums(2, n_s, 3, size(su%t_out)))
    sums = 0.0_dp
    fails = 0

    traj_unit = 0
    if (inp%n_trajectories > 0) then
        open(newunit=traj_unit, file=trim(inp%trajectory_file), status='replace', action='write')
        write(traj_unit, '(a,a)') '# thermal force slab driver, orbit_model: ', trim(inp%orbit_model)
        write(traj_unit, '(a)') '# positions relative to start_x_cm in cm: particle position (full orbit) ' // &
            'or guiding centre (guiding_centre), gradient run'
        write(traj_unit, '(a)') 'particle,t_s,omega_t,x,y,z'
    end if

    !$omp parallel do schedule(dynamic, 64) default(shared) private(ip) reduction(+:sums, fails)
    do ip = 1, inp%n_particles
        call slab_ion(inp, su, ip, traj_unit, sums, fails)
    end do
    !$omp end parallel do

    if (traj_unit /= 0) then
        close(traj_unit)
        print *, 'thermal force: wrote ', trim(inp%trajectory_file)
    end if
    if (any(fails > 0)) print '(a,2i12)', ' thermal force: kernel failures (gradient, control) = ', fails
    call write_slab_output(inp, su, sums, fails)

end subroutine run_slab_driver

subroutine slab_ion(inp, su, ip, traj_unit, sums, fails)

    ! One test ion of the slab driver, gradient and control run in lockstep
    ! with the same random numbers; adds its moments to sums.

    type(thermal_force_input_t), intent(in) :: inp
    type(setup_t), intent(in) :: su
    integer, intent(in) :: ip, traj_unit
    real(dp), intent(inout) :: sums(:,:,:,:)
    integer(int64), intent(inout) :: fails(2)
    type(rng_t) :: rng
    real(dp) :: v(3), x_grad(3), x_ctrl(3), v_grad(3), v_ctrl(3)
    real(dp) :: vpar_grad, vperp_grad, vpar_ctrl, vperp_ctrl
    real(dp) :: xi(3), phase_u, t, t_next
    integer(int64) :: n_step
    integer :: k, n_out, ierr, background_at
    logical :: full_orbit, traced

    n_out = size(su%t_out)
    full_orbit = trim(inp%orbit_model) == 'full_orbit'
    background_at = merge(homma_at_gyro_orbit, homma_at_guiding_centre, &
                          trim(inp%collision_background_at) == 'gyro_orbit')
    traced = ip <= inp%n_trajectories

    call rng_seed(rng, int(inp%seed, int64), int(ip, int64))
    call sample_initial_velocity(inp, su, rng, v)
    v_grad = v
    vpar_grad = 0.0_dp; vperp_grad = 0.0_dp
    ! all ions start with the particle at start_x_cm
    if (full_orbit) then
        x_grad = inp%start_x_cm
    else
        x_grad = inp%start_x_cm - cross(su%b_hat, v)/su%omega_a
        vpar_grad = dot_product(v, su%b_hat)
        vperp_grad = norm2(v - vpar_grad*su%b_hat)
    end if
    x_ctrl = x_grad
    v_ctrl = v_grad
    vpar_ctrl = vpar_grad
    vperp_ctrl = vperp_grad
    t = 0.0_dp
    n_step = 0
    if (traced) call record_trajectory(traj_unit, ip, t, su, x_grad - inp%start_x_cm)

    ! outputs at t = 0 see the initial state
    k = 1
    do while (k <= n_out)
        if (su%t_out(k) > 0.0_dp) exit
        call accumulate(k)
        k = k + 1
    end do

    do while (k <= n_out)
        ! collision at t
        phase_u = rng_uniform(rng)
        call rng_normal3(rng, xi)
        if (full_orbit) then
            call full_orbit_kick(su, su%prof_gradient, su%dt, xi, x_grad, v_grad, ierr)
            if (ierr /= 0) fails(1) = fails(1) + 1
            if (inp%boole_control_run) then
                call full_orbit_kick(su, su%prof_control, su%dt, xi, x_ctrl, v_ctrl, ierr)
                if (ierr /= 0) fails(2) = fails(2) + 1
            end if
        else
            call guiding_centre_kick(su, su%prof_gradient, background_at, inp%boole_collision_shift, su%dt, &
                                     phase_u, xi, x_grad, vpar_grad, vperp_grad, ierr)
            if (ierr /= 0) fails(1) = fails(1) + 1
            if (inp%boole_control_run) then
                call guiding_centre_kick(su, su%prof_control, background_at, inp%boole_collision_shift, su%dt, &
                                         phase_u, xi, x_ctrl, vpar_ctrl, vperp_ctrl, ierr)
                if (ierr /= 0) fails(2) = fails(2) + 1
            end if
            if (traced) call record_trajectory(traj_unit, ip, t, su, x_grad - inp%start_x_cm)
        end if

        ! free motion to the next collision, sampled at the output times
        n_step = n_step + 1
        t_next = n_step*su%dt
        do while (k <= n_out)
            if (su%t_out(k) > t_next) exit
            call free_motion(su%t_out(k) - t)
            t = su%t_out(k)
            call accumulate(k)
            k = k + 1
        end do
        if (k > n_out) exit
        call free_motion(t_next - t)
        t = t_next
    end do

contains

    subroutine free_motion(h)

        real(dp), intent(in) :: h
        integer :: n_sub, i

        if (h <= 0.0_dp) return
        if (full_orbit) then
            if (traced) then
                n_sub = max(1, ceiling(abs(su%omega_a)*h/inp%trajectory_dt_gyro))
                do i = 1, n_sub
                    call helix_push(su, h/n_sub, x_grad, v_grad)
                    call record_trajectory(traj_unit, ip, t + i*h/n_sub, su, x_grad - inp%start_x_cm)
                end do
            else
                call helix_push(su, h, x_grad, v_grad)
            end if
            call helix_push(su, h, x_ctrl, v_ctrl)
        else
            x_grad = x_grad + vpar_grad*h*su%b_hat
            x_ctrl = x_ctrl + vpar_ctrl*h*su%b_hat
        end if

    end subroutine free_motion

    subroutine accumulate(k)

        integer, intent(in) :: k
        real(dp) :: q_grad(n_s), q_ctrl(n_s)

        if (full_orbit) then
            q_grad = [x_grad - inp%start_x_cm, 0.5_dp*su%m_a*dot_product(v_grad, v_grad)]
            q_ctrl = [x_ctrl - inp%start_x_cm, 0.5_dp*su%m_a*dot_product(v_ctrl, v_ctrl)]
        else
            q_grad = [x_grad - inp%start_x_cm, 0.5_dp*su%m_a*(vpar_grad**2 + vperp_grad**2)]
            q_ctrl = [x_ctrl - inp%start_x_cm, 0.5_dp*su%m_a*(vpar_ctrl**2 + vperp_ctrl**2)]
        end if
        sums(1,:,run_gradient,k) = sums(1,:,run_gradient,k) + q_grad
        sums(2,:,run_gradient,k) = sums(2,:,run_gradient,k) + q_grad**2
        sums(1,:,run_control,k) = sums(1,:,run_control,k) + q_ctrl
        sums(2,:,run_control,k) = sums(2,:,run_control,k) + q_ctrl**2
        sums(1,:,run_difference,k) = sums(1,:,run_difference,k) + (q_grad - q_ctrl)
        sums(2,:,run_difference,k) = sums(2,:,run_difference,k) + (q_grad - q_ctrl)**2

    end subroutine accumulate

end subroutine slab_ion

subroutine initialize_gorilla_for_test_ion(inp)

    ! Build GORILLA's mesh and field, then set the test species (mass, charge
    ! and cm_over_e, used by the pusher) and redo the pusher precomputation,
    ! which contains cm_over_e. Requires cylindrical coordinates.

    use constants, only: clight, echarge
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use tetra_physics_mod, only: particle_mass, particle_charge, cm_over_e
    use tetra_physics_poly_precomp_mod, only: make_precomp_poly
    use gorilla_settings_mod, only: coord_system

    type(thermal_force_input_t), intent(in) :: inp

    if (coord_system /= 1) error stop 'thermal force: driver tokamak needs coord_system = 1 (R, phi, Z) in gorilla.inp'
    call initialize_gorilla()
    particle_mass = inp%test_mass_amu*amu_g
    particle_charge = inp%test_charge_number*echarge
    cm_over_e = clight*particle_mass/particle_charge
    call make_precomp_poly()

end subroutine initialize_gorilla_for_test_ion

subroutine tokamak_field(x, ind_tetr, b_field)

    ! Magnetic field vector at x (R, phi, Z) in physical components (R, phi, Z)
    ! from the mesh: covariant h_k = h_k,1 + grad h_k . (x - x_1) per
    ! tetrahedron, physical b = (h_1, h_2/R, h_3) normalised, |B| from bmod_func.

    use tetra_physics_mod, only: tetra_physics
    use supporting_functions_mod, only: bmod_func

    real(dp), intent(in) :: x(3)
    integer, intent(in) :: ind_tetr
    real(dp), intent(out) :: b_field(3)
    real(dp) :: z(3), h(3)

    z = x - tetra_physics(ind_tetr)%x1
    h(1) = tetra_physics(ind_tetr)%h1_1 + sum(tetra_physics(ind_tetr)%gh1*z)
    h(2) = (tetra_physics(ind_tetr)%h2_1 + sum(tetra_physics(ind_tetr)%gh2*z))/x(1)
    h(3) = tetra_physics(ind_tetr)%h3_1 + sum(tetra_physics(ind_tetr)%gh3*z)
    b_field = bmod_func(z, ind_tetr)*h/norm2(h)

end subroutine tokamak_field

subroutine tokamak_start_field(inp, b_field)

    use find_tetra_mod, only: find_tetra

    type(thermal_force_input_t), intent(in) :: inp
    real(dp), intent(out) :: b_field(3)
    real(dp) :: x(3)
    integer :: ind_tetr, iface

    x = inp%start_x_cm          ! find_tetra may map x into its periodic domain
    ! test velocity only for the face disambiguation; v_par must not vanish
    call find_tetra(x, 1.0_dp, 1.0_dp, ind_tetr, iface)
    if (ind_tetr == -1) error stop 'thermal force: start_x_cm lies outside the GORILLA mesh'
    call tokamak_field(x, ind_tetr, b_field)

end subroutine tokamak_start_field

subroutine run_tokamak_driver(inp, su)

    ! As the slab driver: collisions on the fixed grid t_n = n dt, GORILLA
    ! orbit integration in between, sampled at the output times.

    type(thermal_force_input_t), intent(in) :: inp
    type(setup_t), intent(in) :: su
    real(dp), allocatable :: sums(:,:,:,:)
    integer(int64) :: fails(2), n_lost(2)
    integer :: ip

    allocate(sums(2, n_s, 3, size(su%t_out)))
    sums = 0.0_dp
    fails = 0
    n_lost = 0

    !$omp parallel do schedule(dynamic, 16) default(shared) private(ip) reduction(+:sums, fails, n_lost)
    do ip = 1, inp%n_particles
        call tokamak_ion(inp, su, ip, sums, fails, n_lost)
    end do
    !$omp end parallel do

    if (any(fails > 0)) print '(a,2i12)', ' thermal force: kernel failures (gradient, control) = ', fails
    if (any(n_lost > 0)) print '(a,2i12)', ' thermal force: ions lost from the mesh (gradient, control) = ', n_lost
    call write_slab_output(inp, su, sums, fails, n_lost)

end subroutine run_tokamak_driver

subroutine tokamak_ion(inp, su, ip, sums, fails, n_lost)

    ! One test ion on the GORILLA mesh, gradient and control run in lockstep
    ! with the same random numbers. The particle starts at start_x_cm; its
    ! guiding centre at start_x_cm - rho. A lost ion (left the mesh, or the
    ! pusher gave up) keeps its last state and is counted.

    use find_tetra_mod, only: find_tetra
    use anomalous_transport_displacement_mod, only: displace_by_straight_line
    use tetra_grid_settings_mod, only: n_field_periods

    type(thermal_force_input_t), intent(in) :: inp
    type(setup_t), intent(in) :: su
    integer, intent(in) :: ip
    real(dp), intent(inout) :: sums(:,:,:,:)
    integer(int64), intent(inout) :: fails(2), n_lost(2)

    type :: state_t
        real(dp) :: x(3), vpar, vperp
        integer :: ind_tetr, iface, turns
        logical :: lost
    end type state_t

    type(rng_t) :: rng
    type(state_t) :: st(2)
    real(dp) :: v(3), b_hat(3), rho(3), xi(3), phase_u, t, t_next, q(n_s, 2)
    integer(int64) :: n_step
    integer :: k, n_out, r, n_runs, ierr, background_at

    n_out = size(su%t_out)
    n_runs = merge(2, 1, inp%boole_control_run)
    background_at = merge(homma_at_gyro_orbit, homma_at_guiding_centre, &
                          trim(inp%collision_background_at) == 'gyro_orbit')

    call rng_seed(rng, int(inp%seed, int64), int(ip, int64))
    call sample_initial_velocity(inp, su, rng, v)
    b_hat = su%b_hat
    rho = cross(b_hat, v)/su%omega_a
    st(1)%vpar = dot_product(v, b_hat)
    st(1)%vperp = norm2(v - st(1)%vpar*b_hat)
    st(1)%x = inp%start_x_cm - [rho(1), rho(2)/inp%start_x_cm(1), rho(3)]
    st(1)%turns = 0
    st(1)%lost = .false.
    call find_tetra(st(1)%x, st(1)%vpar, st(1)%vperp, st(1)%ind_tetr, st(1)%iface)
    if (st(1)%ind_tetr == -1) st(1)%lost = .true.
    st(2) = st(1)
    do r = 1, 2
        if (st(r)%lost) n_lost(r) = n_lost(r) + 1
    end do

    t = 0.0_dp
    n_step = 0
    k = 1
    do while (k <= n_out)
        if (su%t_out(k) > 0.0_dp) exit
        call accumulate(k)
        k = k + 1
    end do

    do while (k <= n_out)
        ! collision at t, same random numbers for both runs
        phase_u = rng_uniform(rng)
        call rng_normal3(rng, xi)
        call collide(st(1), su%prof_gradient, 1)
        if (n_runs == 2) call collide(st(2), su%prof_control, 2)

        ! GORILLA orbit to the next collision, sampled at the output times
        n_step = n_step + 1
        t_next = n_step*su%dt
        do while (k <= n_out)
            if (su%t_out(k) > t_next) exit
            call push_all(su%t_out(k) - t)
            t = su%t_out(k)
            call accumulate(k)
            k = k + 1
        end do
        if (k > n_out) exit
        call push_all(t_next - t)
        t = t_next
    end do

contains

    subroutine collide(s, prof, r)

        type(state_t), intent(inout) :: s
        type(homma_linear_profile_t), intent(in) :: prof
        integer, intent(in) :: r
        real(dp) :: b_field(3), dx(3), x_local(3)

        if (s%lost) return
        call tokamak_field(s%x, s%ind_tetr, b_field)
        ! local physical frame (R, 0, Z): the profile depends on R and Z only
        x_local = [s%x(1), 0.0_dp, s%x(3)]
        call homma_gc_step(su%units, su%m_a, su%q_a, prof, b_field, background_at, x_local, s%vpar, s%vperp, &
                           su%dt, phase_u, xi, dx, ierr)
        if (ierr /= 0) then
            fails(r) = fails(r) + 1
            return
        end if
        if (.not. inp%boole_collision_shift) return
        call displace_by_straight_line(s%x, s%ind_tetr, s%iface, [dx(1), dx(2)/s%x(1), dx(3)], s%vpar, s%vperp)
        if (s%ind_tetr == -1) then
            s%lost = .true.
            n_lost(r) = n_lost(r) + 1
        end if

    end subroutine collide

    subroutine push_all(h)

        real(dp), intent(in) :: h
        integer :: r

        if (h <= 0.0_dp) return
        do r = 1, n_runs
            call push_orbit(st(r), h, r)
        end do

    end subroutine push_all

    subroutine push_orbit(s, h, r)

        ! GORILLA guiding-centre orbit over h (as orbit_timestep_anomalous_transport)

        use tetra_physics_mod, only: tetra_physics
        use gorilla_settings_mod, only: optional_quantities_type
        use supporting_functions_mod, only: vperp_func
        use utils_orbit_timestep_mod, only: push_within_tetrahedron, initialize_constants_of_motion

        type(state_t), intent(inout) :: s
        real(dp), intent(in) :: h
        integer, intent(in) :: r
        type(optional_quantities_type) :: optional_quantities
        real(dp) :: z_save(3), perpinv, t_remain, t_pass
        integer :: ind_tetr_save, iper_phi
        logical :: boole_t_finished, boole_pusher_failed

        if (s%lost) return
        z_save = s%x - tetra_physics(s%ind_tetr)%x1
        call initialize_constants_of_motion(s%vperp, z_save, s%ind_tetr, perpinv)
        t_remain = h
        do
            ind_tetr_save = s%ind_tetr
            call push_within_tetrahedron(s%ind_tetr, s%iface, s%x, s%vpar, z_save, t_remain, t_pass, &
                                         boole_t_finished, iper_phi, optional_quantities, boole_pusher_failed)
            if (boole_pusher_failed .or. s%ind_tetr == -1) then
                s%lost = .true.
                n_lost(r) = n_lost(r) + 1
                return
            end if
            s%vperp = vperp_func(z_save, perpinv, ind_tetr_save)
            s%turns = s%turns + iper_phi
            t_remain = t_remain - t_pass
            if (boole_t_finished) exit
        end do

    end subroutine push_orbit

    subroutine accumulate(k)

        use constants, only: pi

        integer, intent(in) :: k
        integer :: r
        real(dp) :: phi

        do r = 1, 2
            phi = st(r)%x(2) + st(r)%turns*2.0_dp*pi/n_field_periods
            q(:, r) = [st(r)%x(1) - inp%start_x_cm(1), inp%start_x_cm(1)*(phi - inp%start_x_cm(2)), &
                       st(r)%x(3) - inp%start_x_cm(3), 0.5_dp*su%m_a*(st(r)%vpar**2 + st(r)%vperp**2)]
        end do
        sums(1,:,run_gradient,k) = sums(1,:,run_gradient,k) + q(:,1)
        sums(2,:,run_gradient,k) = sums(2,:,run_gradient,k) + q(:,1)**2
        sums(1,:,run_control,k) = sums(1,:,run_control,k) + q(:,2)
        sums(2,:,run_control,k) = sums(2,:,run_control,k) + q(:,2)**2
        sums(1,:,run_difference,k) = sums(1,:,run_difference,k) + (q(:,1) - q(:,2))
        sums(2,:,run_difference,k) = sums(2,:,run_difference,k) + (q(:,1) - q(:,2))**2

    end subroutine accumulate

end subroutine tokamak_ion

subroutine guiding_centre_kick(su, prof, background_at, shift, h, phase_u, xi, x, vpar, vperp, ierr)

    ! Collision of a guiding centre (homma_gc_step), with its shift if shift.

    type(setup_t), intent(in) :: su
    type(homma_linear_profile_t), intent(in) :: prof
    integer, intent(in) :: background_at
    logical, intent(in) :: shift
    real(dp), intent(in) :: h, phase_u, xi(3)
    real(dp), intent(inout) :: x(3), vpar, vperp
    integer, intent(out) :: ierr
    real(dp) :: dx(3)

    call homma_gc_step(su%units, su%m_a, su%q_a, prof, su%b_field, background_at, x, vpar, vperp, &
                       h, phase_u, xi, dx, ierr)
    if (shift) x = x + dx

end subroutine guiding_centre_kick

subroutine full_orbit_kick(su, prof, h, xi, x, v, ierr)

    ! Collision of a particle at x: Cartesian kick with the background at x.

    type(setup_t), intent(in) :: su
    type(homma_linear_profile_t), intent(in) :: prof
    real(dp), intent(in) :: h, xi(3), x(3)
    real(dp), intent(inout) :: v(3)
    integer, intent(out) :: ierr
    type(homma_local_t) :: loc(1)
    real(dp) :: dv(3)

    call homma_local_at(su%units, su%m_a, su%q_a, prof, x, su%b_field, loc(1), ierr)
    if (ierr /= 0) return
    call homma_kick(loc, v, h, xi, dv, ierr)
    v = v + dv

end subroutine full_orbit_kick

pure subroutine helix_push(su, h, x, v)

    ! Exact motion in the uniform field over h: dv/dt = Omega_a v x b.
    ! v_perp(s) = cos(W s) v_perp - sin(W s) b x v_perp, integrated for x.

    type(setup_t), intent(in) :: su
    real(dp), intent(in) :: h
    real(dp), intent(inout) :: x(3), v(3)
    real(dp) :: v_par(3), v_perp(3), w

    w = su%omega_a
    v_par = dot_product(v, su%b_hat)*su%b_hat
    v_perp = v - v_par
    x = x + v_par*h + sin(w*h)/w*v_perp + (cos(w*h) - 1.0_dp)/w*cross(su%b_hat, v_perp)
    call gyrate(su, -w*h, v)

end subroutine helix_push

subroutine record_trajectory(unit, ip, t, su, x)

    integer, intent(in) :: unit, ip
    real(dp), intent(in) :: t, x(3)
    type(setup_t), intent(in) :: su

    !$omp critical (thermal_force_trajectory)
    write(unit, '(i0,",",4(es24.16,","),es24.16)') ip, t, t*abs(su%omega_a), x
    !$omp end critical (thermal_force_trajectory)

end subroutine record_trajectory

subroutine write_slab_output(inp, su, sums, fails, n_lost)

    ! One row per run (gradient, control, difference) and output time: mean
    ! position relative to start_x_cm (particle for full orbit, guiding centre
    ! otherwise) and kinetic energy, with standard errors (paired for the
    ! difference run).

    type(thermal_force_input_t), intent(in) :: inp
    type(setup_t), intent(in) :: su
    real(dp), intent(in) :: sums(:,:,:,:)
    integer(int64), intent(in) :: fails(2)
    integer(int64), intent(in), optional :: n_lost(2)
    character(len=10), parameter :: run_names(3) = [character(len=10) :: 'gradient', 'control', 'difference']
    real(dp) :: mean(n_s), se(n_s), n
    integer :: unit, k, r, n_runs

    n = real(inp%n_particles, dp)
    n_runs = merge(3, 1, inp%boole_control_run)

    open(newunit=unit, file=trim(inp%output_file), status='replace', action='write')
    write(unit, '(a,a)') '# code: GORILLA_APPLETS thermal_force (i_option = 16), driver ', trim(inp%driver)
    write(unit, '(a,a)') '# orbit_model: ', trim(inp%orbit_model)
    write(unit, '(a,l1)') '# collision_shift: ', inp%boole_collision_shift
    if (trim(inp%driver) == 'tokamak') then
        write(unit, '(a)') '# positions (guiding centre): x = R - R_s, y = R_s (phi - phi_s), z = Z - Z_s; ' // &
            'R_s, phi_s, Z_s = start_x_cm'
        write(unit, '(a,i0,", ",i0)') '# lost_ions_gradient_control: ', n_lost
    else
        write(unit, '(a)') '# positions: relative to start_x_cm; particle (full orbit) or guiding centre'
    end if
    write(unit, '(a,a)') '# collision_background_at: ', trim(inp%collision_background_at)
    write(unit, '(a,a)') '# heat_flux_model: ', trim(inp%heat_flux_model)
    write(unit, '(a,3es24.16)') '# grad_temperature_ev_cm: ', inp%grad_temperature_ev_cm
    write(unit, '(a,es24.16,a,es24.16)') '# test_mass_amu: ', inp%test_mass_amu, ', test_charge_number: ', &
        inp%test_charge_number
    write(unit, '(a,es24.16)') '# tau_s_s: ', su%tau_s
    write(unit, '(a,es24.16)') '# omega_a_1_s: ', su%omega_a
    write(unit, '(a,es24.16,a,es24.16)') '# dt_coll_s: ', su%dt, ', dt_coll_over_tau_s: ', inp%dt_coll_over_tau_s
    write(unit, '(a,i0,a,i0)') '# n_particles: ', inp%n_particles, ', seed: ', inp%seed
    write(unit, '(a,i0,", ",i0)') '# kernel_failures_gradient_control: ', fails
    write(unit, '(a)') '# units: CGS; positions in cm relative to start_x_cm, energy in erg'
    write(unit, '(a)') 'run,t_s,t_over_tau_s,omega_t,mean_x,se_x,mean_y,se_y,mean_z,se_z,mean_energy,se_energy'

    do k = 1, size(su%t_out)
        do r = 1, n_runs
            mean = sums(1,:,r,k)/n
            se = sqrt(max(sums(2,:,r,k)/n - mean**2, 0.0_dp)/(n - 1.0_dp))
            write(unit, '(a,",",11(es24.16,:,","))') trim(run_names(r)), su%t_out(k), su%t_out(k)/su%tau_s, &
                su%t_out(k)*abs(su%omega_a), mean(s_x), se(s_x), mean(s_y), se(s_y), mean(s_z), se(s_z), &
                mean(s_energy), se(s_energy)
        end do
    end do
    close(unit)
    print *, 'thermal force: wrote ', trim(inp%output_file)

end subroutine write_slab_output

subroutine sample_initial_velocity(inp, su, rng, v)

    use constants, only: ev2erg

    type(thermal_force_input_t), intent(in) :: inp
    type(setup_t), intent(in) :: su
    type(rng_t), intent(inout) :: rng
    real(dp), intent(out) :: v(3)
    real(dp) :: xi(3)

    call rng_normal3(rng, xi)
    select case (trim(inp%test_distribution))
    case ('maxwellian')
        v = sqrt(inp%test_temperature_ev*ev2erg/su%m_a)*xi
    case ('drifting_bimaxwellian')
        v = sqrt(inp%test_temperature_par_ev*ev2erg/su%m_a)*xi(1)*su%b_hat &
          + sqrt(inp%test_temperature_perp_ev*ev2erg/su%m_a)*(xi(2)*su%e1 + xi(3)*su%e2)
    end select
    v = v + inp%test_drift_cm_s

end subroutine sample_initial_velocity

pure subroutine gyrate(su, angle, v)

    ! Rotate the perpendicular velocity about b: exact solution of
    ! dv/dt = Omega_a v x b over a time angle/(-Omega_a).

    type(setup_t), intent(in) :: su
    real(dp), intent(in) :: angle
    real(dp), intent(inout) :: v(3)
    real(dp) :: v_par(3), v_perp(3)

    v_par = dot_product(v, su%b_hat)*su%b_hat
    v_perp = v - v_par
    v = v_par + cos(angle)*v_perp + sin(angle)*cross(su%b_hat, v_perp)

end subroutine gyrate

pure function moments(su, v) result(q)

    type(setup_t), intent(in) :: su
    real(dp), intent(in) :: v(3)
    real(dp) :: q(n_q), v_par, v2

    v_par = dot_product(v, su%b_hat)
    v2 = dot_product(v, v)
    q(q_vx:q_vz) = v
    q(q_vpar) = v_par
    q(q_vpar2) = v_par**2
    q(q_vperp2) = v2 - v_par**2
    q(q_v4) = v2**2

end function moments

subroutine write_velocity_output(inp, su, sums, fails, ev2erg)

    ! One row per run (gradient, control, difference) and output time.
    ! Per test ion (density 1): p_par = m_a (<v_par^2> - <v_par>^2),
    ! p_perp = m_a <v_perp^2>/2, v4 = <|v|^4>. Standard errors are those of
    ! the per-ion sample means; for the difference run they are paired.

    type(thermal_force_input_t), intent(in) :: inp
    type(setup_t), intent(in) :: su
    real(dp), intent(in) :: sums(:,:,:,:), ev2erg
    integer(int64), intent(in) :: fails(2)
    character(len=10), parameter :: run_names(3) = [character(len=10) :: 'gradient', 'control', 'difference']
    real(dp) :: mean(n_q, 3), se(n_q, 3), p_par(3), p_perp(3), n
    integer :: unit, k, r, n_runs

    n = real(inp%n_particles, dp)
    n_runs = merge(3, 1, inp%boole_control_run)

    open(newunit=unit, file=trim(inp%output_file), status='replace', action='write')
    write(unit, '(a)') '# code: GORILLA_APPLETS thermal_force (i_option = 16), driver velocity'
    write(unit, '(a)') '# operator: Homma & Hatayama 2014, Euler-Maruyama in Cartesian velocity, exact gyration between kicks'
    write(unit, '(a,a)') '# heat_flux_model: ', trim(inp%heat_flux_model)
    write(unit, '(a,3es24.16)') '# grad_temperature_ev_cm: ', inp%grad_temperature_ev_cm
    write(unit, '(a,es24.16,a,es24.16)') '# test_mass_amu: ', inp%test_mass_amu, ', test_charge_number: ', &
        inp%test_charge_number
    write(unit, '(a,es24.16)') '# tau_s_s: ', su%tau_s
    write(unit, '(a,es24.16)') '# omega_a_1_s: ', su%omega_a
    write(unit, '(a,es24.16,a,es24.16)') '# dt_coll_s: ', su%dt, ', dt_coll_over_tau_s: ', inp%dt_coll_over_tau_s
    write(unit, '(a,i0,a,i0)') '# n_particles: ', inp%n_particles, ', seed: ', inp%seed
    write(unit, '(a,l1)') '# gyration: ', inp%boole_gyration
    write(unit, '(a,i0,", ",i0)') '# kernel_failures_gradient_control: ', fails
    write(unit, '(a,es24.16)') '# ev2erg: ', ev2erg
    write(unit, '(a)') '# units: CGS; velocities cm/s, p_par and p_perp per test ion in erg, v4 in cm^4/s^4'
    write(unit, '(a)') 'run,t_s,t_over_tau_s,omega_t,mean_vx,se_vx,mean_vy,se_vy,mean_vz,se_vz,' // &
        'u_par,se_u_par,p_par,se_p_par,p_perp,se_p_perp,v4,se_v4'

    do k = 1, size(su%t_out)
        do r = 1, 3
            mean(:,r) = sums(1,:,r,k)/n
            se(:,r) = sqrt(max(sums(2,:,r,k)/n - mean(:,r)**2, 0.0_dp)/(n - 1.0_dp))
        end do
        do r = 1, 2
            p_par(r) = su%m_a*(mean(q_vpar2,r) - mean(q_vpar,r)**2)
            p_perp(r) = 0.5_dp*su%m_a*mean(q_vperp2,r)
        end do
        p_par(3) = p_par(1) - p_par(2)
        p_perp(3) = p_perp(1) - p_perp(2)
        do r = 1, n_runs
            write(unit, '(a,",",17(es24.16,:,","))') trim(run_names(r)), su%t_out(k), su%t_out(k)/su%tau_s, &
                su%t_out(k)*abs(su%omega_a), mean(q_vx,r), se(q_vx,r), mean(q_vy,r), se(q_vy,r), &
                mean(q_vz,r), se(q_vz,r), mean(q_vpar,r), se(q_vpar,r), p_par(r), su%m_a*se(q_vpar2,r), &
                p_perp(r), 0.5_dp*su%m_a*se(q_vperp2,r), mean(q_v4,r), se(q_v4,r)
        end do
    end do
    close(unit)
    print *, 'thermal force: wrote ', trim(inp%output_file)

end subroutine write_velocity_output

pure function cross(x, y) result(z)

    real(dp), intent(in) :: x(3), y(3)
    real(dp) :: z(3)

    z = [x(2)*y(3) - x(3)*y(2), x(3)*y(1) - x(1)*y(3), x(1)*y(2) - x(2)*y(1)]

end function cross

end module thermal_force_mod
