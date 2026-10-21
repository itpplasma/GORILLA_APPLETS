program test_deltaf_kernels
    ! Independent oracles for the FortSym-generated delta-f kernels:
    ! centred finite differences of the generated values reproduce the generated
    ! derivatives, and the circular chart agrees with GORILLA's analytic circular
    ! field (field_analytic_circ_mod, grid_kind = 5).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use deltaf_phase_step_kernel, only: evaluate_deltaf_phase_step
    use deltaf_phase_derivative_kernel, only: evaluate_deltaf_phase_derivative
    use deltaf_complex_harmonic_kernel, only: evaluate_deltaf_complex_harmonic
    use deltaf_harmonic_kernel, only: evaluate_deltaf_harmonic
    use deltaf_circular_chart_kernel, only: evaluate_deltaf_circular_chart
    use deltaf_weight_kernel, only: evaluate_deltaf_weight
    use field_analytic_circ_mod, only: set_field_analytic_circ, field_analytic_circ, &
        psi_pol_analytic_circ
    implicit none

    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    ! Deuterium in Gaussian units.
    real(dp), parameter :: mass = 3.3435837724e-24_dp, charge = 4.80320471e-10_dp
    real(dp), parameter :: c_light = 2.99792458e10_dp
    ! Circular tokamak (GORILLA field_analytic_circ_mod parameters).
    real(dp), parameter :: r0 = 165.0_dp, a = 50.0_dp, b0 = 2.0e4_dp
    real(dp), parameter :: q0 = 1.1_dp, q1 = 2.4_dp
    integer :: failures

    failures = 0
    call check_harmonic()
    call check_complex_harmonic()
    call check_phase_recurrence()
    call check_circular_against_gorilla()
    call check_weight()
    if (failures /= 0) then
        write (*, "(a,i0,a)") "test_deltaf_kernels: ", failures, " failures"
        error stop 1
    end if
    write (*, "(a)") "test_deltaf_kernels: all checks passed"

contains

    subroutine check_phase_recurrence()
        real(dp) :: sn, cs, sn_next, cs_next, th, ph, angle, dph, dth, refph, refth
        integer :: j, m

        do j = 0, 16
            th = j*2*pi/16
            ph = 0.132_dp + j*0.7_dp
            do m = -320, 320
                if (mod(m + 320, 32) == 0) then
                    angle = m*th - 3*ph
                    sn = sin(angle)
                    cs = cos(angle)
                end if
                call evaluate_deltaf_phase_derivative(real(m, dp), 3.0_dp, &
                    1.3_dp, -0.7_dp, sn, cs, dph, dth)
                call evaluate_deltaf_complex_harmonic(real(m, dp), 3.0_dp, &
                    1.3_dp, -0.7_dp, th, ph, refph, refth)
                call check_close("phase recurrence phi", dph, refph, 1.0e-10_dp)
                call check_close("phase recurrence theta", dth, refth, 1.0e-10_dp)
                call evaluate_deltaf_phase_step(sn, cs, sin(th), cos(th), sn_next, cs_next)
                sn = sn_next
                cs = cs_next
            end do
        end do
    end subroutine check_phase_recurrence

    subroutine check_complex_harmonic()
        real(dp) :: m, n, bc, bs, th, ph, alpha, dph, dth
        integer :: i

        do i = -320, 320
            m = real(i, dp)
            n = 3.0_dp
            bc = 1.3_dp
            bs = -0.7_dp
            th = 0.37_dp
            ph = 1.41_dp
            alpha = m*th - n*ph
            call evaluate_deltaf_complex_harmonic(m, n, bc, bs, th, ph, dph, dth)
            call check_close("complex dphi", dph, &
                n*(bc*sin(alpha) - bs*cos(alpha)), 1.0e-12_dp)
            call check_close("complex dtheta", dth, &
                m*(-bc*sin(alpha) + bs*cos(alpha)), 1.0e-12_dp)
        end do
    end subroutine check_complex_harmonic

    subroutine check_harmonic()
        real(dp), parameter :: m_pol = 2.0_dp, n_tor = 3.0_dp, h = 1.0e-6_dp
        real(dp) :: v(8), vp(8), vm(8), s0, th, ph
        real(dp) :: bmn, dbmn, phimn, dphimn

        s0 = 0.31_dp
        th = 0.7_dp
        ph = 1.9_dp
        call profiles(s0, bmn, dbmn, phimn, dphimn)
        call evaluate_deltaf_harmonic(m_pol, n_tor, bmn, dbmn, phimn, dphimn, th, &
            ph, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8))
        call check_close("delta B value", v(1), bmn*cos(m_pol*th - n_tor*ph), 1.0e-14_dp)
        call profiles(s0 + h, bmn, dbmn, phimn, dphimn)
        call evaluate_deltaf_harmonic(m_pol, n_tor, bmn, dbmn, phimn, dphimn, th, &
            ph, vp(1), vp(2), vp(3), vp(4), vp(5), vp(6), vp(7), vp(8))
        call profiles(s0 - h, bmn, dbmn, phimn, dphimn)
        call evaluate_deltaf_harmonic(m_pol, n_tor, bmn, dbmn, phimn, dphimn, th, &
            ph, vm(1), vm(2), vm(3), vm(4), vm(5), vm(6), vm(7), vm(8))
        call check_close("d deltaB/ds", v(2), (vp(1) - vm(1))/(2*h), 1.0e-7_dp)
        call check_close("d deltaPhi/ds", v(6), (vp(5) - vm(5))/(2*h), 1.0e-7_dp)
        call profiles(s0, bmn, dbmn, phimn, dphimn)
        call evaluate_deltaf_harmonic(m_pol, n_tor, bmn, dbmn, phimn, dphimn, &
            th + h, ph, vp(1), vp(2), vp(3), vp(4), vp(5), vp(6), vp(7), vp(8))
        call evaluate_deltaf_harmonic(m_pol, n_tor, bmn, dbmn, phimn, dphimn, &
            th - h, ph, vm(1), vm(2), vm(3), vm(4), vm(5), vm(6), vm(7), vm(8))
        call check_close("d deltaB/dtheta", v(3), (vp(1) - vm(1))/(2*h), 1.0e-7_dp)
        call check_close("d deltaPhi/dtheta", v(7), (vp(5) - vm(5))/(2*h), 1.0e-7_dp)
        call evaluate_deltaf_harmonic(m_pol, n_tor, bmn, dbmn, phimn, dphimn, th, &
            ph + h, vp(1), vp(2), vp(3), vp(4), vp(5), vp(6), vp(7), vp(8))
        call evaluate_deltaf_harmonic(m_pol, n_tor, bmn, dbmn, phimn, dphimn, th, &
            ph - h, vm(1), vm(2), vm(3), vm(4), vm(5), vm(6), vm(7), vm(8))
        call check_close("d deltaB/dphi", v(4), (vp(1) - vm(1))/(2*h), 1.0e-7_dp)
        call check_close("d deltaPhi/dphi", v(8), (vp(5) - vm(5))/(2*h), 1.0e-7_dp)
    end subroutine check_harmonic

    !> Test harmonic profiles: b_mn = 2 G sqrt(s), Phi_mn = 3 statV s (1 - s).
    subroutine profiles(s, bmn, dbmn, phimn, dphimn)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: bmn, dbmn, phimn, dphimn

        bmn = 2.0_dp*sqrt(s)
        dbmn = 1.0_dp/sqrt(s)
        phimn = 3.0_dp*s*(1.0_dp - s)
        dphimn = 3.0_dp*(1.0_dp - 2.0_dp*s)
    end subroutine profiles

    subroutine check_circular_against_gorilla()
        real(dp), parameter :: h = 1.0e-4_dp
        real(dp) :: radius, th, out(11), outp(11), outm(11)
        real(dp) :: rr, zz, br, bp, bz, dbrdr, dbrdp, dbrdz, dbpdr, dbpdp, dbpdz
        real(dp) :: dbzdr, dbzdp, dbzdz, bmod_g, dbmod_dr_g
        integer :: i, j

        call set_field_analytic_circ(r0, a, b0, q0, q1, 2)
        do i = 1, 4
            do j = 1, 5
                radius = 0.2_dp*a*i
                th = 0.37_dp + 1.21_dp*j
                call circ(radius, th, out)
                rr = r0 + radius*cos(th)
                zz = radius*sin(th)
                call field_analytic_circ(rr, 0.4_dp, zz, br, bp, bz, dbrdr, dbrdp, &
                    dbrdz, dbpdr, dbpdp, dbpdz, dbzdr, dbzdp, dbzdz)
                bmod_g = sqrt(br**2 + bp**2 + bz**2)
                dbmod_dr_g = (br*(dbrdr*cos(th) + dbrdz*sin(th)) &
                    + bp*(dbpdr*cos(th) + dbpdz*sin(th)) &
                    + bz*(dbzdr*cos(th) + dbzdz*sin(th)))/bmod_g
                call check_close("circular |B| vs GORILLA", out(1), bmod_g, 1.0e-13_dp)
                call check_close("circular d|B|/dr vs GORILLA", out(2), dbmod_dr_g, &
                    1.0e-11_dp)
                call check_close("circular psi_pol vs GORILLA psi_pol", out(4), &
                    psi_pol_analytic_circ(radius), 1.0e-13_dp)
                call check_close("circular h_phi vs GORILLA R B_phi/B", out(8), &
                    rr*bp/bmod_g, 1.0e-13_dp)
                call check_close("circular R_major", out(11), rr, 1.0e-15_dp)
                call circ(radius + h, th, outp)
                call circ(radius - h, th, outm)
                call check_close("circular dpsi/dr FD", out(5), (outp(4) - outm(4))/(2*h), &
                    1.0e-8_dp)
                call check_close("circular ds/dr FD", out(10), (outp(9) - outm(9))/(2*h), &
                    1.0e-8_dp)
                call check_close("circular dpsi/ds", out(6), out(5)/out(10), 1.0e-13_dp)
                call circ(radius, th + h, outp)
                call circ(radius, th - h, outm)
                call check_close("circular d|B|/dtheta FD", out(3), &
                    (outp(1) - outm(1))/(2*h), 1.0e-7_dp)
            end do
        end do
    end subroutine check_circular_against_gorilla

    subroutine circ(rad, ang, o)
        real(dp), intent(in) :: rad, ang
        real(dp), intent(out) :: o(11)

        call evaluate_deltaf_circular_chart(r0, a, b0, q0, q1, rad, ang, o(1), &
            o(2), o(3), o(4), o(5), o(6), o(7), o(8), o(9), o(10), o(11))
    end subroutine circ

    !> Canonical Maxwellian with n = n0 exp(-psi/L_n), T = T0 (1 + psi/L_T),
    !> Phi0 = E0 psi. d ln f0/dp_phi is compared with a centred difference of
    !> the generated ln f0 in p_phi at fixed H0.
    subroutine check_weight()
        real(dp), parameter :: mu = 2.0e-12_dp, vpar = -2.4e7_dp, bmod = 2.1e4_dp
        real(dp), parameter :: db1 = 3.0_dp, dphi1 = 0.02_dp
        real(dp) :: p0, hp, h0, out(5), outp(5), outm(5)

        h0 = 0.5_dp*mass*vpar**2 + mu*bmod + charge*0.4_dp
        p0 = charge*(-2.0e6_dp)/c_light + mass*vpar*0.9_dp*165.0_dp
        hp = 1.0e-6_dp*abs(p0)
        call weight_at(p0, h0, out)
        call weight_at(p0 + hp, h0, outp)
        call weight_at(p0 - hp, h0, outm)
        call check_close("d ln f0/dp_phi FD", out(2), (outp(1) - outm(1))/(2*hp), 1.0e-6_dp)
        call check_close("dH1/dphi = (mu + m v_par^2/B) dB1/dphi + e dPhi1/dphi", out(3), &
            (mu + mass*vpar**2/bmod)*db1 + charge*dphi1, 1.0e-14_dp)
        call check_close("dw/dt = dlnf0/dp_phi dH1/dphi", out(4), out(2)*out(3), 1.0e-14_dp)
        call check_close("torque sample = -dH1/dphi", out(5), -out(3), 1.0e-15_dp)
    end subroutine check_weight

    subroutine weight_at(p, h0, o)
        real(dp), intent(in) :: p, h0
        real(dp), intent(out) :: o(5)
        real(dp), parameter :: mu = 2.0e-12_dp, vpar = -2.4e7_dp, bmod = 2.1e4_dp
        real(dp), parameter :: db1 = 3.0_dp, dphi1 = 0.02_dp
        real(dp), parameter :: n0 = 5.0e13_dp, ln = 3.0e7_dp, t0 = 1.6e-9_dp
        real(dp), parameter :: lt = 2.0e7_dp, e0 = -1.0e-7_dp
        real(dp) :: psi, dens, temp

        psi = c_light*p/charge
        dens = n0*exp(-psi/ln)
        temp = t0*(1.0_dp + psi/lt)
        call evaluate_deltaf_weight(mass, charge, c_light, mu, vpar, bmod, h0, &
            dens, -dens/ln, temp, t0/lt, e0*psi, e0, db1, dphi1, o(1), o(2), &
            o(3), o(4), o(5))
    end subroutine weight_at

    subroutine check_close(label, got, want, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, rtol

        call close_abs(label, got, want, rtol*max(abs(want), tiny(1.0_dp)))
    end subroutine check_close

    subroutine close_abs(label, got, want, atol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, atol

        if (abs(got - want) <= atol) return
        failures = failures + 1
        write (*, "(a,3es24.15)") "FAIL "//label//" got/want/tol:", got, want, atol
    end subroutine close_abs

end program test_deltaf_kernels
