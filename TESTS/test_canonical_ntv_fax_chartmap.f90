! Independent circular-coordinate oracle for the chartmap f_AX mirror source.
program test_canonical_ntv_fax_chartmap
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use canonical_ntv_mod, only: fax_surface_field, fax_surface_b, &
        ntv_harmonic_dphi, ntv_m_pol, ntv_n_tor, ntv_bmn, ntv_amplitude, &
        ntv_bmn_relative
    use tetra_grid_settings_mod, only: netcdf_filename, load_tetra_grid_inp
    use boozer_chartmap_mod, only: load_boozer_chartmap, evaluate_boozer_chartmap, &
        chartmap_torflux, chartmap_iota
    use deltaf_fax_chartmap_kernel, only: evaluate_deltaf_fax_chartmap
    use deltaf_fax_weight_kernel, only: evaluate_deltaf_fax_weight
    use deltaf_fax_coordinates_kernel, only: evaluate_deltaf_fax_coordinates
    implicit none
    real(dp), parameter :: r0 = 165.0_dp, a = 69.991288_dp, b0 = -17991.7336_dp
    real(dp), parameter :: mass = 3.3452e-24_dp
    real(dp), parameter :: mu = 1.0e-13_dp, vpar = 1.0e7_dp, kappa = -0.0023_dp
    character(len=512) :: writer
    real(dp) :: x(3), s, r, q, th, tg, rmaj, bexact, bmin, bmax, blo, bhi
    real(dp) :: at, ap, az, hs, ht, hp, b, jac, ct, cp, iota, dbt, dbp, dpot
    real(dp) :: dw, old, ref, gradpar, gbth, lam, vt, lam0, vt0, dchi, dchi0
    real(dp) :: relative_db1
    real(dp) :: err_b, err_vt, err_m0, err_m3, err_deriv, phase, denom
    integer :: i, m, rc, u

    call get_command_argument(1, writer)
    call execute_command_line(trim(writer)//' fax_axi.nc 0 10 3 2 8193 128 2', &
        exitstat=rc)
    if (rc /= 0) error stop 'chartmap writer failed'
    open(newunit=u, file='tetra_grid.inp', status='replace')
    write(u, '(a)') "&tetra_grid_nml grid_kind=6, netcdf_filename='fax_axi.nc' /"
    close(u)
    call load_tetra_grid_inp()
    call load_boozer_chartmap(trim(netcdf_filename))
    call fax_surface_field()
    err_b = 0.0_dp
    err_vt = 0.0_dp
    err_m0 = 0.0_dp
    err_m3 = 0.0_dp
    err_deriv = 0.0_dp
    ntv_amplitude = 1.0_dp
    ntv_bmn = 10.0_dp
    ntv_n_tor = 2
    do i = 1, 20
        s = 0.07_dp + 0.041_dp*real(i, dp)
        th = (-1.0_dp)**i*(0.3_dp + 0.11_dp*real(i, dp))
        x = [s, th, 0.231_dp*real(i, dp)]
        r = sqrt(r0**2 - (r0 - s*(r0 - sqrt(r0**2 - a**2)))**2)
        q = 1.0_dp + 2.0_dp*(r/a)**2
        tg = 2.0_dp*atan2(sqrt(1.0_dp + r/r0)*sin(th/2.0_dp), &
            sqrt(1.0_dp - r/r0)*cos(th/2.0_dp))
        rmaj = r0 + r*cos(tg)
        bexact = abs(b0)/rmaj*sqrt(r0**2 + (r/q)**2)
        blo = abs(b0)/(r0 + r)*sqrt(r0**2 + (r/q)**2)
        bhi = abs(b0)/(r0 - r)*sqrt(r0**2 + (r/q)**2)
        call fax_surface_b(x, bmin, bmax)
        err_b = max(err_b, abs(bmin/blo - 1.0_dp), abs(bmax/bhi - 1.0_dp))
        call evaluate_boozer_chartmap(s, th, x(3), az, at, ap, hs, ht, hp, b, jac)
        ct = chartmap_iota(s)*chartmap_torflux
        cp = chartmap_torflux
        ! Exact circular b^phi and b^theta (not the generated normalization).
        iota = sqrt(1.0_dp - (r/r0)**2)/q
        denom = b0*r0/(rmaj**2*bexact)
        gbth = abs(b0)*sqrt(r0**2 + (r/q)**2)*r*sin(tg)/rmaj**2 &
            * (1.0_dp + r/r0*cos(tg))/sqrt(1.0_dp - (r/r0)**2)
        call evaluate_deltaf_fax_chartmap(mu, ht, hp, ct, cp, gbth, 0.0_dp, &
            kappa, 0.0_dp, 0.0_dp, 1.0_dp, b, 0.0_dp, dw, gradpar)
        if (gradpar*th >= 0.0_dp) &
            error stop 'chartmap branch does not follow b.grad B'
        call evaluate_deltaf_fax_coordinates(mass, mu, vpar, &
            min(max(b, bmin), bmax), bmin, bmax, sign(1.0_dp, gradpar), lam, vt)
        call evaluate_deltaf_fax_coordinates(mass, mu, vpar, bexact, blo, bhi, &
            -sign(1.0_dp, th), lam0, vt0)
        err_vt = max(err_vt, abs(vt - vt0))
        dchi = 0.3_dp + 0.1_dp*vt + 0.2_dp*lam
        dchi0 = 0.3_dp + 0.1_dp*vt0 + 0.2_dp*lam0
        do m = 0, 3, 3
            ntv_m_pol = m
            call ntv_harmonic_dphi(s, th, x(3), dbp, dpot, b, dbt)
            phase = real(m, dp)*th - 2.0_dp*x(3)
            err_deriv = max(err_deriv, abs(dbp - 20.0_dp*sin(phase)), &
                abs(dbt + real(m, dp)*10.0_dp*sin(phase)))
            call evaluate_deltaf_fax_chartmap(mu, ht, hp, ct, cp, gbth, 0.0_dp, &
                kappa, dbt, dbp, dchi, b, 0.0_dp, dw, gradpar)
            if (m == 0) then
                call evaluate_deltaf_fax_weight(mu, b0*r0/bexact, rmaj, &
                    kappa, dbp, dchi0, old)
                err_m0 = max(err_m0, abs(dw/old - 1.0_dp))
            else
                ! Analytic helical parallel derivative: (n - m*iota) b^phi b sin.
                ref = kappa*mu*(b0*r0/bexact)*denom &
                    * (2.0_dp - real(m, dp)*iota)*10.0_dp*sin(phase)*dchi0
                err_m3 = max(err_m3, abs(dw/ref - 1.0_dp))
            end if
        end do
    end do
    print '(a, es12.4)', 'surface extrema max relative error: ', err_b
    print '(a, es12.4)', 'vartheta max absolute error [rad]: ', err_vt
    print '(a, es12.4)', 'm=0 chartmap / cylindrical f_AX max relative error: ', err_m0
    print '(a, es12.4)', 'm=3 chartmap / analytic mirror max relative error: ', err_m3
    print '(a, es12.4)', 'harmonic angular derivatives max absolute error: ', err_deriv
    if (err_b > 2.0e-3_dp .or. err_vt > 1.0e-2_dp .or. &
        err_m0 > 5.0e-3_dp .or. err_m3 > 1.0e-2_dp .or. &
        err_deriv > 1.0e-12_dp) error stop 'circular f_AX oracle failed'

    ! Sign reversal and a resonant nonzero-m perturbation (iota=2/3, m=3,n=2).
    call evaluate_deltaf_fax_chartmap(mu, -2.0_dp, -3.0_dp, 2.0_dp, 3.0_dp, &
        1.0_dp, 2.0_dp, kappa, -3.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, dw, gradpar)
    if (abs(dw) > 0.0_dp .or. gradpar >= 0.0_dp) error stop 'resonance/branch failed'
    call evaluate_deltaf_fax_chartmap(mu, 2.0_dp, 3.0_dp, 2.0_dp, 3.0_dp, &
        1.0_dp, 2.0_dp, kappa, -3.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, dw, gradpar)
    if (abs(dw) > 0.0_dp .or. gradpar <= 0.0_dp) error stop 'reversed branch failed'
    ! Relative amplitudes also differentiate the local B factor.
    ntv_bmn_relative = .true.
    call ntv_harmonic_dphi(s, th, x(3), dbp, dpot, b, dbt, relative_db1)
    if (abs(relative_db1/(10.0_dp*b*cos(phase)) - 1.0_dp) > 1.0e-13_dp) &
        error stop 'relative harmonic value failed'
    call evaluate_deltaf_fax_chartmap(1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 3.0_dp, &
        5.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 10.0_dp, 2.0_dp, &
        dw, gradpar)
    if (abs(dw - 42.0_dp/130.0_dp) > 1.0e-14_dp) &
        error stop 'relative amplitude mirror derivative failed'
    print '(a)', 'PASS: circular f_AX chartmap and nonzero-m analytic oracles'
end program test_canonical_ntv_fax_chartmap
