!Follow one guiding-centre orbit on a Boozer chartmap mesh (grid_kind 6) with
!the vertex potential set by set_flux_mesh_potential and print the rate of the
!field-line label alpha = theta - iota*phi, the radial cell of the start point
!and the final position (hex, for bitwise comparisons).
!Usage: test_flux_mesh_potential_driver.x <mode> <s0> <energy_eV> <t_step> <n_steps> <p0>
!  mode none:    GORILLA default potential (eps_Phi*A_theta)
!  mode zero:    set_flux_mesh_potential with Phi0 = 0
!  mode profile: set_flux_mesh_potential with Phi0(s) = p0*(s - 0.8*s**2) [statV]
module flux_mesh_potential_test_profile

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    real(dp) :: p0 = 0.0_dp

contains

    function phi_profile(s) result(phi)
        real(dp), intent(in) :: s
        real(dp) :: phi
        phi = p0*(s - 0.8_dp*s**2)
    end function phi_profile

    function phi_zero(s) result(phi)
        real(dp), intent(in) :: s
        real(dp) :: phi
        phi = 0.0_dp*s
    end function phi_zero

end module flux_mesh_potential_test_profile

program test_flux_mesh_potential_driver

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use constants, only: pi, ev2erg, clight
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp, n_field_periods
    use tetra_grid_mod, only: tetra_grid, verts_sthetaphi
    use tetra_physics_mod, only: particle_mass
    use orbit_timestep_gorilla_mod, only: initialize_gorilla, orbit_timestep_gorilla
    use boozer_chartmap_mod, only: chartmap_iota, chartmap_torflux
    use flux_mesh_potential_mod, only: set_flux_mesh_potential
    use flux_mesh_potential_test_profile, only: p0, phi_profile, phi_zero

    implicit none

    character(len=64) :: arg, mode
    real(dp) :: s0, energy_ev, t_step, vmod, vpar, vperp, x(3), theta_old, phi_old
    real(dp) :: theta_unwrapped, phi_unwrapped, iota, smin, smax, s_lo, s_hi
    real(dp) :: t, alpha, st, sa, stt, sta, rate
    integer :: n_steps, k, ind_tetr, iface
    logical :: boole_initialized

    if (command_argument_count() /= 6) error stop &
        'usage: <mode> <s0> <energy_eV> <t_step> <n_steps> <p0>'
    call get_command_argument(1, mode)
    call get_command_argument(2, arg)
    read(arg, *) s0
    call get_command_argument(3, arg)
    read(arg, *) energy_ev
    call get_command_argument(4, arg)
    read(arg, *) t_step
    call get_command_argument(5, arg)
    read(arg, *) n_steps
    call get_command_argument(6, arg)
    read(arg, *) p0

    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla(ipert_in=0)
    select case (trim(mode))
    case ('none')
    case ('zero')
        call set_flux_mesh_potential(phi_zero, 0)
    case ('profile')
        call set_flux_mesh_potential(phi_profile, 0)
    case default
        error stop 'mode must be none, zero or profile'
    end select

    !Deeply trapped start on the outboard midplane (|B| minimum). A small
    !nonzero v_par: the pusher divides by v_par for its reference time.
    vmod = sqrt(2.0_dp*energy_ev*ev2erg/particle_mass)
    vpar = 1.0e-3_dp*vmod
    vperp = sqrt(vmod**2 - vpar**2)
    x = [s0, 0.0_dp, 0.1_dp]
    iota = chartmap_iota(s0)

    theta_old = x(2)
    phi_old = x(3)
    theta_unwrapped = x(2)
    phi_unwrapped = x(3)
    smin = s0
    smax = s0
    st = 0.0_dp
    sa = 0.0_dp
    stt = 0.0_dp
    sta = 0.0_dp
    s_lo = -1.0_dp
    s_hi = -1.0_dp
    boole_initialized = .false.
    ind_tetr = -1
    iface = -1
    do k = 1, n_steps
        call orbit_timestep_gorilla(x, vpar, vperp, t_step, boole_initialized, ind_tetr, iface)
        if (ind_tetr == -1) error stop 'test orbit left the mesh'
        if (k == 1) then
            s_lo = minval(verts_sthetaphi(1, tetra_grid(ind_tetr)%ind_knot))
            s_hi = maxval(verts_sthetaphi(1, tetra_grid(ind_tetr)%ind_knot))
        end if
        theta_unwrapped = theta_unwrapped + wrap(x(2) - theta_old, 2.0_dp*pi)
        phi_unwrapped = phi_unwrapped &
                        + wrap(x(3) - phi_old, 2.0_dp*pi/real(n_field_periods, dp))
        theta_old = x(2)
        phi_old = x(3)
        smin = min(smin, x(1))
        smax = max(smax, x(1))
        t = real(k, dp)*t_step
        alpha = theta_unwrapped - iota*phi_unwrapped
        st = st + t
        sa = sa + alpha
        stt = stt + t*t
        sta = sta + t*alpha
    end do
    rate = (real(n_steps, dp)*sta - st*sa)/(real(n_steps, dp)*stt - st*st)

    print '(a,es24.16)', 'RESULT alpha_rate ', rate
    print '(a,es24.16)', 'RESULT s0 ', s0
    print '(a,es24.16)', 'RESULT s_lo ', s_lo
    print '(a,es24.16)', 'RESULT s_hi ', s_hi
    print '(a,es24.16)', 'RESULT s_min ', smin
    print '(a,es24.16)', 'RESULT s_max ', smax
    print '(a,es24.16)', 'RESULT iota ', iota
    print '(a,es24.16)', 'RESULT torflux ', chartmap_torflux
    print '(a,es24.16)', 'RESULT clight ', clight
    print '(a,es24.16)', 'RESULT theta_turns ', (theta_unwrapped - 0.0_dp)/(2.0_dp*pi)
    print '(a,3z17)', 'FINAL ', transfer(x, 0_int64, 3)

contains

    pure function wrap(delta, period) result(wrapped)
        real(dp), intent(in) :: delta, period
        real(dp) :: wrapped
        wrapped = modulo(delta + 0.5_dp*period, period) - 0.5_dp*period
    end function wrap

end program test_flux_mesh_potential_driver
