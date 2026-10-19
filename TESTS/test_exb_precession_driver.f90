!Follow one guiding-centre orbit on the flux-tube start position of the
!mono-energetic transport applet and print the rate of change of the
!field-line label alpha = theta - iota*phi (Boozer angles of a chartmap).
!Parallel motion leaves alpha unchanged, so the difference of two runs with
!v_E /= 0 and v_E = 0 isolates the ExB rotation.
!Usage: test_exb_precession_driver.x <pitch> <t_step> <n_steps>
program test_exb_precession_driver

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use constants, only: pi, ev2erg
    use gorilla_settings_mod, only: load_gorilla_inp, eps_Phi
    use tetra_grid_settings_mod, only: load_tetra_grid_inp, n_field_periods
    use gorilla_applets_settings_mod, only: load_gorilla_applets_inp
    use mono_energetic_transp_coef_settings_mod, only: load_mono_energetic_transp_coef_inp, &
                                                       energy_eV
    use gorilla_applets_sub_mod, only: initialize_mono_energetic_transp_coef
    use fluxtv_mod, only: pos_fluxtv_mat
    use tetra_physics_mod, only: particle_mass
    use orbit_timestep_gorilla_mod, only: orbit_timestep_gorilla
    use boozer_chartmap_mod, only: chartmap_iota

    implicit none

    character(len=64) :: arg
    real(dp) :: pitch, t_step, vmod, vpar, vperp, x(3), theta_old, phi_old
    real(dp) :: theta_unwrapped, phi_unwrapped, iota, s0, smin, smax
    real(dp) :: t, alpha, st, sa, stt, sta, rate
    integer :: n_steps, k, ind_tetr, iface
    logical :: boole_initialized

    if (command_argument_count() /= 3) error stop 'usage: <pitch> <t_step> <n_steps>'
    call get_command_argument(1, arg)
    read(arg, *) pitch
    call get_command_argument(2, arg)
    read(arg, *) t_step
    call get_command_argument(3, arg)
    read(arg, *) n_steps

    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call load_gorilla_applets_inp()
    call load_mono_energetic_transp_coef_inp()
    call initialize_mono_energetic_transp_coef()

    vmod = sqrt(2.0_dp*energy_eV*ev2erg/particle_mass)
    vpar = pitch*vmod
    vperp = sqrt(1.0_dp - pitch**2)*vmod
    x = pos_fluxtv_mat(1, 1:3)
    s0 = x(1)
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
    boole_initialized = .false.
    ind_tetr = -1
    iface = -1
    do k = 1, n_steps
        call orbit_timestep_gorilla(x, vpar, vperp, t_step, boole_initialized, ind_tetr, iface)
        if (ind_tetr == -1) error stop 'test orbit left the mesh'
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
    print '(a,es24.16)', 'RESULT eps_Phi ', eps_Phi
    print '(a,es24.16)', 'RESULT vmod ', vmod
    print '(a,es24.16)', 'RESULT s0 ', s0
    print '(a,es24.16)', 'RESULT s_excursion ', smax - smin
    print '(a,es24.16)', 'RESULT phi_turns ', (phi_unwrapped - pos_fluxtv_mat(1, 3))/(2.0_dp*pi)
    print '(a,es24.16)', 'RESULT theta_turns ', (theta_unwrapped - pos_fluxtv_mat(1, 2))/(2.0_dp*pi)
    print '(a,es24.16)', 'RESULT iota ', iota

contains

    pure function wrap(delta, period) result(wrapped)
        real(dp), intent(in) :: delta, period
        real(dp) :: wrapped
        wrapped = modulo(delta + 0.5_dp*period, period) - 0.5_dp*period
    end function wrap

end program test_exb_precession_driver
