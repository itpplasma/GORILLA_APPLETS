! Native orbit oracle: deterministic pitch ensemble, optional Lorentz kicks.
! Input chart/mesh is supplied in the working directory. This tests the same
! GORILLA pusher used by the estimator, independently of its weight equation.
program deltaf_3d_orbit_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_timestep_gorilla_mod, only: initialize_gorilla, orbit_timestep_gorilla
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp, n_field_periods
    use tetra_physics_mod, only: particle_mass
    use canonical_ntv_mod, only: ntv_lorentz_proposal
    use constants, only: ev2erg, pi
    implicit none
    integer :: n, steps, i, j, ind, iface, flips, last_sign, current_sign, seedn, u
    integer, allocatable :: seed(:)
    real(dp) :: x(3), xp(3), vpar, vperp, speed, energy, dt, nu, trace, lam
    real(dp) :: r(4), xi, lam_new, log_q, theta_delta, zeta_delta, sd, t0, t1
    real(dp) :: s_start, theta_start, phi_start, pitch, poloidal, toroidal
    real(dp) :: phi_min, phi_max
    logical :: initialized, valid
    character(len=128) :: arg
    if (command_argument_count() /= 10) &
        error stop 'usage: probe n steps dt nu E_eV s theta phi pitch seed'
    call get_command_argument(1, arg); read (arg, *) n
    call get_command_argument(2, arg); read (arg, *) steps
    call get_command_argument(3, arg); read (arg, *) dt
    call get_command_argument(4, arg); read (arg, *) nu
    call get_command_argument(5, arg); read (arg, *) energy
    call get_command_argument(6, arg); read (arg, *) s_start
    call get_command_argument(7, arg); read (arg, *) theta_start
    call get_command_argument(8, arg); read (arg, *) phi_start
    call get_command_argument(9, arg); read (arg, *) pitch
    call get_command_argument(10, arg); read (arg, *) seedn
    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla()
    call random_seed(size=i)
    allocate (seed(i))
    seed = [(seedn + 7919*j, j=1, i)]
    call random_seed(put=seed)
    speed = sqrt(2*energy*ev2erg/particle_mass)
    open (newunit=u, file='probe.dat', status='replace')
    write (u, '(a)') &
        '# n lost s_start s_end theta_travel phi_travel flips trace_s phi_span'
    call cpu_time(t0)
    do j = 1, n
        call random_number(r)
        x = [s_start, theta_start, phi_start]
        if (n > 1) x(2:3) = 2*pi*r(1:2)
        lam = pitch
        if (n > 1) lam = 2*r(3) - 1
        vpar = speed*lam
        vperp = speed*sqrt(1 - lam**2)
        initialized = .false.
        ind = 0
        iface = 0
        flips = 0
        poloidal = 0
        toroidal = 0
        phi_min = 0
        phi_max = 0
        last_sign = int(sign(1.0_dp, vpar))
        trace = 0
        do i = 1, steps
            xp = x
            call orbit_timestep_gorilla(x, vpar, vperp, dt, initialized, ind, iface)
            if (ind == -1) exit
            theta_delta = modulo(x(2) - xp(2) + pi, 2*pi) - pi
            zeta_delta = modulo(x(3) - xp(3) + pi/real(n_field_periods, dp), &
                          2*pi/real(n_field_periods, dp)) - pi/real(n_field_periods, dp)
            poloidal = poloidal + theta_delta
            toroidal = toroidal + zeta_delta
            phi_min = min(phi_min, toroidal)
            phi_max = max(phi_max, toroidal)
            current_sign = int(sign(1.0_dp, vpar))
            if (current_sign /= last_sign) flips = flips + 1
            last_sign = current_sign
            trace = trace + dt
            if (nu > 0) then
                call random_number(r)
                xi = sqrt(-2*log(max(r(1), tiny(1.0_dp))))*cos(2*pi*r(2))
                sd = sqrt(vpar**2 + vperp**2)
                lam = vpar/sd
                call ntv_lorentz_proposal(lam, nu*dt, xi, lam_new, log_q, valid)
                if (valid) then
                    vpar = sd*lam_new
                    vperp = sd*sqrt(1 - lam_new**2)
                end if
            end if
        end do
        write (u, '(2i8,4es24.15,i8,2es24.15)') j, merge(1, 0, ind == -1), &
            s_start, x(1), poloidal, toroidal, flips, trace, phi_max - phi_min
    end do
    call cpu_time(t1)
    close (u)
    print '(a,es24.15)', 'PROBE seconds per marker ', (t1 - t0)/real(n, dp)
end program deltaf_3d_orbit_probe
