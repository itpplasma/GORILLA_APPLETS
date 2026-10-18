program test_deltaf_kernels
    ! Independent behavioural oracles for the FortSym-generated delta-f kernels
    ! in SRC/generated (generator: gen/deltaf_kernels).
    !   1. Centred finite differences of the generated values reproduce the
    !      generated derivatives (harmonic, circular chart, flux label,
    !      canonical momentum, canonical-Maxwellian d ln f0/dp_phi).
    !   2. The circular chart agrees with GORILLA's analytic circular field
    !      (field_analytic_circ_mod, grid_kind = 5): |B|, psif, covariant h_phi
    !      and d|B|/dr from GORILLA's cylindrical field derivatives.
    !   3. A manufactured Boozer chartmap (grid_kind = 6) read back through
    !      GORILLA's boozer_chartmap_mod: interpolated Bmod reproduces the
    !      generated delta B, a finite difference of interpolated Bmod in zeta
    !      reproduces the generated dH1/dphi, and GORILLA's A_phi, h_phi give
    !      the generated canonical momentum.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use deltaf_harmonic_kernel, only: evaluate_deltaf_harmonic
    use deltaf_circular_chart_kernel, only: evaluate_deltaf_circular_chart
    use deltaf_circular_flux_label_kernel, only: evaluate_deltaf_circular_flux_label
    use deltaf_canonical_momentum_kernel, only: evaluate_deltaf_canonical_momentum
    use deltaf_weight_kernel, only: evaluate_deltaf_weight
    use field_analytic_circ_mod, only: set_field_analytic_circ, field_analytic_circ
    use field_eq_mod, only: psif
    use boozer_chartmap_mod, only: load_boozer_chartmap, evaluate_boozer_chartmap
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
    call check_circular_against_gorilla()
    call check_flux_label()
    call check_canonical_momentum()
    call check_weight()
    call check_boozer_chartmap()
    if (failures /= 0) then
        write (*, "(a,i0,a)") "test_deltaf_kernels: ", failures, " failures"
        error stop 1
    end if
    write (*, "(a)") "test_deltaf_kernels: all checks passed"

contains

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

        call set_field_analytic_circ(r0, a, b0, q0, q1)
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
                call check_close("circular psi_pol vs GORILLA psif", out(4), psif, 1.0e-13_dp)
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

    subroutine check_flux_label()
        real(dp) :: out(11), s_star, ds_dpsi, sp, sm, dummy, hpsi

        call evaluate_deltaf_circular_chart(r0, a, b0, q0, q1, 23.0_dp, 0.3_dp, &
            out(1), out(2), out(3), out(4), out(5), out(6), out(7), out(8), out(9), &
            out(10), out(11))
        call evaluate_deltaf_circular_flux_label(a, b0, q0, q1, out(4), s_star, &
            ds_dpsi)
        call check_close("flux label inverts psi_pol(s)", s_star, out(9), 1.0e-13_dp)
        call check_close("ds*/dpsi* = 1/(dpsi/ds)", ds_dpsi, 1.0_dp/out(6), 1.0e-12_dp)
        hpsi = 1.0e-6_dp*abs(out(4))
        call evaluate_deltaf_circular_flux_label(a, b0, q0, q1, out(4) + hpsi, sp, &
            dummy)
        call evaluate_deltaf_circular_flux_label(a, b0, q0, q1, out(4) - hpsi, sm, &
            dummy)
        call check_close("ds*/dpsi* FD", ds_dpsi, (sp - sm)/(2*hpsi), 1.0e-7_dp)
    end subroutine check_flux_label

    subroutine check_canonical_momentum()
        real(dp) :: p, psi_star, dp_dv, pp, pm, dummy1, dummy2, vpar, h

        vpar = 3.1e7_dp
        h = 1.0e2_dp
        call evaluate_deltaf_canonical_momentum(mass, charge, c_light, vpar, &
            -1.3e6_dp, 4.1e6_dp, 2.2e4_dp, p, psi_star, dp_dv)
        call check_close("p_phi = e A_phi/c + m v_par B_phi/B", p, &
            charge*(-1.3e6_dp)/c_light + mass*vpar*4.1e6_dp/2.2e4_dp, 1.0e-14_dp)
        call check_close("psi* = c p_phi/e", psi_star, c_light*p/charge, 1.0e-14_dp)
        call evaluate_deltaf_canonical_momentum(mass, charge, c_light, vpar + h, &
            -1.3e6_dp, 4.1e6_dp, 2.2e4_dp, pp, dummy1, dummy2)
        call evaluate_deltaf_canonical_momentum(mass, charge, c_light, vpar - h, &
            -1.3e6_dp, 4.1e6_dp, 2.2e4_dp, pm, dummy1, dummy2)
        call check_close("dp_phi/dv_par FD", dp_dv, (pp - pm)/(2*h), 1.0e-6_dp)
    end subroutine check_canonical_momentum

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

    subroutine check_boozer_chartmap()
        character(*), parameter :: fixture = "deltaf_kernels_boozer_chartmap.nc"
        integer, parameter :: nfp = 3
        real(dp), parameter :: m_pol = 2.0_dp, n_tor = 3.0_dp
        real(dp), parameter :: mu = 2.0e-12_dp, vpar = 1.7e7_dp
        real(dp) :: s, th, ze, a_s, a_t, a_p, h_s, h_t, h_p, bmod_g, sqrtg
        real(dp) :: v(8), bmn, dbmn, phimn, dphimn, b_axi, dz, bp, bm, dummy(8)
        real(dp) :: p_gen, psi_star, dp_dv, dh1(5), dw
        integer :: i

        call write_boozer_fixture(fixture, nfp)
        call load_boozer_chartmap(fixture)
        do i = 1, 6
            s = 0.12_dp + 0.13_dp*i
            th = 0.41_dp + 0.97_dp*i
            ze = 0.23_dp + 0.61_dp*i
            call evaluate_boozer_chartmap(s, th, ze, a_s, a_t, a_p, h_s, h_t, h_p, &
                bmod_g, sqrtg)
            call booz_profiles(s, b_axi, bmn, dbmn)
            phimn = 0.0_dp
            dphimn = 0.0_dp
            call evaluate_deltaf_harmonic(m_pol, n_tor, bmn, dbmn, phimn, dphimn, &
                th, ze, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8))
            call close_abs("Boozer delta B vs GORILLA interpolated Bmod", v(1), &
                bmod_g - b_axi, 5.0e-3_dp*bmn)
            ! Interpolated Bmod in zeta, centred over two grid cells.
            dz = 2.0_dp*(2.0_dp*pi/nfp)/64.0_dp
            call evaluate_boozer_chartmap(s, th, ze + dz, dummy(1), dummy(2), &
                dummy(3), dummy(4), dummy(5), dummy(6), bp, dummy(7))
            call evaluate_boozer_chartmap(s, th, ze - dz, dummy(1), dummy(2), &
                dummy(3), dummy(4), dummy(5), dummy(6), bm, dummy(7))
            call evaluate_deltaf_weight(mass, charge, c_light, mu, vpar, bmod_g, &
                1.0e-8_dp, 1.0e13_dp, 0.0_dp, 1.6e-9_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                v(4), v(8), dh1(1), dh1(2), dh1(3), dw, dh1(5))
            call close_abs("Boozer dH1/dphi vs GORILLA FD in zeta", dh1(3), &
                (mu + mass*vpar**2/bmod_g)*(bp - bm)/(2*dz), &
                2.0e-2_dp*(mu + mass*vpar**2/bmod_g)*n_tor*bmn)
            call evaluate_deltaf_canonical_momentum(mass, charge, c_light, vpar, &
                a_p, h_p*bmod_g, bmod_g, p_gen, psi_star, dp_dv)
            call check_close("Boozer p_phi vs GORILLA A_phi, h_phi", p_gen, &
                charge*a_p/c_light + mass*vpar*h_p, 1.0e-13_dp)
            call check_close("Boozer GORILLA A_phi(s) is the fixture profile", a_p, &
                -0.4_dp*2.0e8_dp*s, 1.0e-12_dp)
        end do
        open (newunit=i, file=fixture, status="old")
        close (i, status="delete")
    end subroutine check_boozer_chartmap

    !> Axisymmetric |B| and harmonic amplitude, linear in rho = sqrt(s) so the
    !> chartmap's radial interpolation is exact and only the harmonic's angular
    !> interpolation error remains (b_axi is theta-independent for that reason).
    subroutine booz_profiles(s, b_axi, bmn, dbmn)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: b_axi, bmn, dbmn

        b_axi = 2.0e4_dp*(1.0_dp - 0.1_dp*sqrt(s))
        bmn = 40.0_dp*sqrt(s)
        dbmn = 20.0_dp/sqrt(s)
    end subroutine booz_profiles

    subroutine write_boozer_fixture(filename, nfp)
        character(*), intent(in) :: filename
        integer, intent(in) :: nfp
        integer, parameter :: nr = 33, nt = 128, nz = 64
        integer :: ncid, dims(4), vars(12), ir, it, iz, dim1(1), dim3(3)
        real(dp) :: rho(nr), s(nr), theta(nt), zeta(nz), radius, b_axi, bmn, dbmn
        real(dp), allocatable :: x(:, :, :), y(:, :, :), zz(:, :, :), bmod(:, :, :)
        real(dp) :: aphi(nr), btheta(nr), bphi(nr)

        allocate (x(nr, nt, nz), y(nr, nt, nz), zz(nr, nt, nz), bmod(nr, nt, nz))

        do ir = 1, nr
            rho(ir) = real(ir - 1, dp)/real(nr - 1, dp)
        end do
        s = rho**2
        do it = 1, nt
            theta(it) = 2.0_dp*pi*real(it - 1, dp)/real(nt, dp)
        end do
        do iz = 1, nz
            zeta(iz) = (2.0_dp*pi/nfp)*real(iz - 1, dp)/real(nz, dp)
        end do
        aphi = -0.4_dp*2.0e8_dp*s
        btheta = 1.0e5_dp*rho
        bphi = 3.3e6_dp - 1.0e4_dp*rho
        do iz = 1, nz
            do it = 1, nt
                do ir = 1, nr
                    radius = 165.0_dp + 50.0_dp*rho(ir)*cos(theta(it))
                    x(ir, it, iz) = radius*cos(zeta(iz))
                    y(ir, it, iz) = radius*sin(zeta(iz))
                    zz(ir, it, iz) = 50.0_dp*rho(ir)*sin(theta(it))
                    call booz_profiles(s(ir), b_axi, bmn, dbmn)
                    bmod(ir, it, iz) = b_axi + 40.0_dp*rho(ir) &
                        *cos(2.0_dp*theta(it) - 3.0_dp*zeta(iz))
                end do
            end do
        end do

        call nc(nf90_create(filename, nf90_clobber, ncid))
        call nc(nf90_def_dim(ncid, "rho", nr, dims(1)))
        call nc(nf90_def_dim(ncid, "s", nr, dims(2)))
        call nc(nf90_def_dim(ncid, "theta", nt, dims(3)))
        call nc(nf90_def_dim(ncid, "zeta", nz, dims(4)))
        dim1 = dims(1)
        call nc(nf90_def_var(ncid, "rho", nf90_double, dim1, vars(1)))
        call nc(nf90_def_var(ncid, "B_theta", nf90_double, dim1, vars(6)))
        call nc(nf90_def_var(ncid, "B_phi", nf90_double, dim1, vars(7)))
        dim1 = dims(2)
        call nc(nf90_def_var(ncid, "s", nf90_double, dim1, vars(2)))
        call nc(nf90_def_var(ncid, "A_phi", nf90_double, dim1, vars(5)))
        dim1 = dims(3)
        call nc(nf90_def_var(ncid, "theta", nf90_double, dim1, vars(3)))
        dim1 = dims(4)
        call nc(nf90_def_var(ncid, "zeta", nf90_double, dim1, vars(4)))
        dim3 = dims([1, 3, 4])
        call nc(nf90_def_var(ncid, "x", nf90_double, dim3, vars(8)))
        call nc(nf90_def_var(ncid, "y", nf90_double, dim3, vars(9)))
        call nc(nf90_def_var(ncid, "z", nf90_double, dim3, vars(10)))
        call nc(nf90_def_var(ncid, "Bmod", nf90_double, dim3, vars(11)))
        call nc(nf90_def_var(ncid, "num_field_periods", nf90_int, vars(12)))
        call nc(nf90_put_att(ncid, nf90_global, "torflux", 2.0e8_dp))
        call nc(nf90_put_att(ncid, nf90_global, "boozer_field", 1))
        call nc(nf90_enddef(ncid))
        call nc(nf90_put_var(ncid, vars(1), rho))
        call nc(nf90_put_var(ncid, vars(2), s))
        call nc(nf90_put_var(ncid, vars(3), theta))
        call nc(nf90_put_var(ncid, vars(4), zeta))
        call nc(nf90_put_var(ncid, vars(5), aphi))
        call nc(nf90_put_var(ncid, vars(6), btheta))
        call nc(nf90_put_var(ncid, vars(7), bphi))
        call nc(nf90_put_var(ncid, vars(8), x))
        call nc(nf90_put_var(ncid, vars(9), y))
        call nc(nf90_put_var(ncid, vars(10), zz))
        call nc(nf90_put_var(ncid, vars(11), bmod))
        call nc(nf90_put_var(ncid, vars(12), nfp))
        call nc(nf90_close(ncid))
    end subroutine write_boozer_fixture

    subroutine nc(status)
        integer, intent(in) :: status

        if (status /= nf90_noerr) error stop nf90_strerror(status)
    end subroutine nc

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
