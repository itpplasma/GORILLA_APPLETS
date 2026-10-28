!> Unit checks of the Boozer chartmap route of the canonical NTV estimator on
!> the chartmap of GORILLA's analytic circular field (write_circular_boozer_
!> chartmap.x, path in argument 1):
!>   1. the chartmap reproduces the analytic field (A_phi, B_phi, |B|, R, Z, to
!>      the accuracy of GORILLA's linear node interpolation) and
!>      its Jacobian formula matches the geometric Jacobian (Boozer chart);
!>   2. Fourier projection recovers the cosine and sine harmonics exactly;
!>   3. the flux table (psi_pol = A_phi(s), psi_tor = torflux s, q) agrees
!>      with GORILLA's analytic flux_functions table: psi -> s inverse;
!>   4. the J0 window volume equals the analytic torus-annulus volume and the
!>      zeta average B0 removes the perturbation;
!>   5. the wall polygon test.
program test_canonical_ntv_chartmap

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use chartmap_ntv_mod, only: chartmap_raw_t, read_chartmap_raw, chartmap_harmonic, &
        write_chartmap_flux_functions, chartmap_mesh_setup, chartmap_point, &
        chartmap_window_volume
    use boozer_chartmap_mod, only: load_boozer_chartmap, evaluate_boozer_chartmap, &
        evaluate_chartmap_geometry
    use field_analytic_circ_mod, only: set_field_analytic_circ, &
        psi_pol_analytic_circ, q_analytic_circ
    use wall_contour_mod, only: load_wall_contour, wall_contains, clear_wall_contour

    implicit none

    real(dp), parameter :: r0 = 165.0_dp, a = 69.991288_dp, b0 = -17991.7336_dp
    real(dp), parameter :: pi = acos(-1.0_dp), bmn = 12.0_dp, b2 = 5.0_dp
    character(len=*), parameter :: fixture = 'circular_boozer_chartmap.nc'
    character(len=512) :: writer
    integer :: nfail = 0
    real(dp) :: s_edge

    call get_command_argument(1, writer)
    call execute_command_line(trim(writer)//' '//fixture// &
        ' 1.0 12.0 2 3 33 64 48 -1 3 5.0', exitstat=nfail)
    if (nfail /= 0) error stop 'fixture writer failed'
    call set_field_analytic_circ(r0, a, b0, 1.0_dp, 2.0_dp, 2)
    s_edge = r0 - sqrt(r0**2 - a**2)

    call check_field()
    call check_harmonics()
    call check_flux_table()
    call check_volume()
    call check_wall()
    if (nfail > 0) then
        print '(a, i0, a)', 'FAIL: ', nfail, ' checks'
        error stop 1
    end if
    print '(a)', 'PASS: canonical NTV chartmap unit checks'

contains

    real(dp) function r_of_s(s)
        real(dp), intent(in) :: s
        r_of_s = sqrt(max(0.0_dp, r0**2 - (r0 - s*s_edge)**2))
    end function r_of_s

    real(dp) function theta_geo(ths, r)
        real(dp), intent(in) :: ths, r
        theta_geo = 2.0_dp*atan2(sqrt(1.0_dp + r/r0)*sin(0.5_dp*ths), &
                                 sqrt(1.0_dp - r/r0)*cos(0.5_dp*ths))
    end function theta_geo

    subroutine check_field()
        real(dp) :: s, th, ze, a_s, a_t, a_p, h_s, h_t, h_p, b, sqrtg, r, tg, rr, pc, zz
        real(dp) :: drds, dzds, bexact, jgeo, dthg
        integer :: i

        call load_boozer_chartmap(fixture)
        do i = 1, 7
            s = 0.1_dp + 0.12_dp*i
            th = 0.37_dp + 0.83_dp*i
            ze = 0.21_dp + 0.5_dp*i
            r = r_of_s(s)
            tg = theta_geo(th, r)
            call evaluate_boozer_chartmap(s, th, ze, a_s, a_t, a_p, h_s, h_t, h_p, b, sqrtg)
            call evaluate_chartmap_geometry(s, th, ze, rr, pc, zz, drds, dzds)
            bexact = abs(b0)/(r0 + r*cos(tg))*sqrt(r0**2 + (r/q_analytic_circ(r))**2) &
                + bmn*cos(2.0_dp*th - 3.0_dp*ze) + b2*sin(-th - 3.0_dp*ze)
            call close_rel('A_phi', a_p, psi_pol_analytic_circ(r), 1.0e-3_dp)
            call close_rel('B_phi = h_phi |B|', h_p*b, b0*r0, 1.0e-12_dp)
            call close_rel('|B|', b, bexact, 2.0e-3_dp)
            call close_abs('R', rr, r0 + r*cos(tg), 0.6_dp)
            call close_abs('Z', zz, r*sin(tg), 0.6_dp)
            ! Geometric Jacobian of (s, theta*, phi): -r R dr/ds dtheta/dtheta*.
            dthg = (1.0_dp + r/r0*cos(tg))/sqrt(1.0_dp - (r/r0)**2)
            jgeo = -r*(r0 + r*cos(tg))*(r0*sqrt(1.0_dp - (r/r0)**2)*s_edge/r)*dthg
            call close_rel('Boozer Jacobian vs geometry', sqrtg, jgeo, 5.0e-3_dp)
        end do
    end subroutine check_field

    subroutine check_harmonics()
        type(chartmap_raw_t) :: cm
        real(dp), allocatable :: bc(:), bs(:)

        call read_chartmap_raw(fixture, cm, with_geometry=.false.)
        allocate(bc(cm%nrho), bs(cm%nrho))
        call chartmap_harmonic(cm, 2, 3, bc, bs)
        call close_abs('bc(2,3)', maxval(abs(bc - bmn)), 0.0_dp, 1.0e-9_dp)
        call close_abs('bs(2,3)', maxval(abs(bs)), 0.0_dp, 1.0e-9_dp)
        call chartmap_harmonic(cm, -1, 3, bc, bs)
        call close_abs('bc(-1,3)', maxval(abs(bc)), 0.0_dp, 1.0e-9_dp)
        call close_abs('bs(-1,3)', maxval(abs(bs - b2)), 0.0_dp, 1.0e-9_dp)
        call chartmap_harmonic(cm, 1, 3, bc, bs)
        call close_abs('(1,3) absent', maxval(abs(bc)) + maxval(abs(bs)), 0.0_dp, 1.0e-9_dp)
        call chartmap_harmonic(cm, 2, 2, bc, bs)
        call close_abs('(2,2) absent', maxval(abs(bc)) + maxval(abs(bs)), 0.0_dp, 1.0e-9_dp)
    end subroutine check_harmonics

    !> Flux table of the chartmap against GORILLA's analytic table: on the
    !> analytic rows, s from the chartmap A_phi(s) inverse (piecewise linear,
    !> as GORILLA interpolates A_phi) equals psi_tor/psi_tor_edge.
    subroutine check_flux_table()
        type(chartmap_raw_t) :: cm
        real(dp) :: row(5), psi(200), stab(200), qtab(200), reff(200), s_inv, r, err
        integer :: u, ios, n, i, k
        character(len=256) :: line

        call read_chartmap_raw(fixture, cm)
        call write_chartmap_flux_functions(cm, 'chartmap_flux_functions.dat')
        open(newunit=u, file='chartmap_flux_functions.dat', status='old')
        n = 0
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '#') > 0) cycle
            read(line, *) row
            n = n + 1
            psi(n) = row(4)
            stab(n) = row(5)
            qtab(n) = row(3)
            reff(n) = row(2)
        end do
        close(u)
        if (n /= cm%ns) call fail('flux table row count')
        stab(:n) = stab(:n)/stab(n)
        call close_abs('table s = psi_tor/psi_tor_edge', maxval(abs(stab(:n) - cm%s)), &
            0.0_dp, 1.0e-14_dp)
        call close_abs('table psi_pol = A_phi', maxval(abs(psi(:n) - cm%aphi)), 0.0_dp, &
            1.0e-6_dp*maxval(abs(cm%aphi)))
        err = 0.0_dp
        do i = 2, n - 1
            r = r_of_s(cm%s(i))
            err = max(err, abs(qtab(i)*sqrt(1.0_dp - (r/r0)**2)/q_analytic_circ(r) - 1.0_dp))
            err = max(err, abs(reff(i) - r)/a)
        end do
        call close_abs('table q (winding) and r_eff', err, 0.0_dp, 5.0e-3_dp)
        err = 0.0_dp
        do i = 1, 99
            r = a*real(i, dp)/100.0_dp
            k = 1
            do while (k < n - 1 .and. (psi(k + 1) - psi_pol_analytic_circ(r)) &
                      *(psi(n) - psi(1)) < 0.0_dp)
                k = k + 1
            end do
            s_inv = cm%s(k) + (cm%s(k + 1) - cm%s(k))*(psi_pol_analytic_circ(r) - psi(k)) &
                /(psi(k + 1) - psi(k))
            err = max(err, abs(s_inv - (r0 - sqrt(r0**2 - r**2))/s_edge))
        end do
        call close_abs('psi -> s inverse vs analytic s(psi)', err, 0.0_dp, 2.0e-3_dp)
    end subroutine check_flux_table

    subroutine check_volume()
        real(dp) :: vol, jmax, r1, r2, b, b00, j0, bphi, aphi, r, tg, bexact
        integer :: i

        call chartmap_mesh_setup(fixture, 'chartmap_flux_functions.dat')
        call chartmap_window_volume(0.3_dp, 0.8_dp, 2, vol, jmax)
        r1 = r_of_s(0.3_dp)
        r2 = r_of_s(0.8_dp)
        call close_rel('J0 window volume (one of 2 periods)', vol, &
            pi*pi*r0*(r2**2 - r1**2), 2.0e-3_dp)
        do i = 1, 5
            call chartmap_point(0.2_dp*i - 0.05_dp, 0.7_dp*i, 0.9_dp*i, b, b00, j0, &
                bphi, aphi)
            r = r_of_s(0.2_dp*i - 0.05_dp)
            tg = theta_geo(0.7_dp*i, r)
            bexact = abs(b0)/(r0 + r*cos(tg))*sqrt(r0**2 + (r/q_analytic_circ(r))**2)
            call close_rel('B0 = zeta average', b00, bexact, 2.0e-3_dp)
            call close_abs('B - B0 = perturbation', b - b00, &
                bmn*cos(1.4_dp*i - 2.7_dp*i) + b2*sin(-0.7_dp*i - 2.7_dp*i), 0.5_dp)
        end do
    end subroutine check_volume

    subroutine check_wall()
        integer :: u

        open(newunit=u, file='wall.dat', status='replace')
        write(u, '(a)') '# L-shaped (concave) limiter, R Z [cm]'
        write(u, '(a)') '100 -50'
        write(u, '(a)') '200 -50'
        write(u, '(a)') '200 0'
        write(u, '(a)') '150 0'
        write(u, '(a)') '150 50'
        write(u, '(a)') '100 50'
        close(u)
        call load_wall_contour('wall.dat')
        if (.not. wall_contains(120.0_dp, 20.0_dp)) call fail('wall: inside arm')
        if (.not. wall_contains(180.0_dp, -20.0_dp)) call fail('wall: inside foot')
        if (wall_contains(180.0_dp, 20.0_dp)) call fail('wall: concave notch')
        if (wall_contains(90.0_dp, 0.0_dp)) call fail('wall: outside')
        call clear_wall_contour()
        if (.not. wall_contains(1.0e4_dp, 0.0_dp)) call fail('wall: inactive')
    end subroutine check_wall

    subroutine close_rel(label, got, want, rtol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: got, want, rtol
        if (abs(got - want) > rtol*abs(want)) then
            print '(a, a, 2es23.15)', 'mismatch ', label, got, want
            nfail = nfail + 1
        end if
    end subroutine close_rel

    subroutine close_abs(label, got, want, atol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: got, want, atol
        if (abs(got - want) > atol) then
            print '(a, a, 2es23.15)', 'mismatch ', label, got, want
            nfail = nfail + 1
        end if
    end subroutine close_abs

    subroutine fail(label)
        character(len=*), intent(in) :: label
        print '(2a)', 'failed: ', label
        nfail = nfail + 1
    end subroutine fail

end program test_canonical_ntv_chartmap
