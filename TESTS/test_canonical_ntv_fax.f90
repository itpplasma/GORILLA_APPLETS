!> Collisional background drive of the canonical NTV estimator (ntv_fax_file):
!> (a) the generated weight source for h = G (chi = lambda, so
!>     D = d chi/d eta_hat = -B/(2 Bmax lambda)) equals the mu part
!>     (d ln f0/dp_phi) mu (h_phi/R)**2 d dB/dphi of the canonical source, i.e.
!>     the mirror force on the canonical shift (generator identity F3);
!> (b) the generated coordinates: eta_hat = 1 on the trapped-passing boundary
!>     lambda**2 = 1 - B/Bmax, vartheta = 0, branch pi/2, branch pi at
!>     B = Bmin, (Bmin + Bmax)/2, Bmax;
!> (c) the table reader and bilinear lookup reproduce a bilinear D exactly and
!>     clamp outside the table.
program test_canonical_ntv_fax

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use deltaf_fax_coordinates_kernel, only: evaluate_deltaf_fax_coordinates
    use deltaf_fax_weight_kernel, only: evaluate_deltaf_fax_weight
    use canonical_ntv_mod, only: ntv_fax_load_table, ntv_fax_dchi_deta

    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), mass = 3.3452e-24_dp
    real(dp), parameter :: bmin = 1.74e4_dp, bmax = 2.25e4_dp
    integer :: k, nfail
    real(dp) :: tau, b, vpar, mu, v, lam, eta, vt, hphi, rmaj, kappa, db1, dw, ref
    real(dp) :: err, err_max, d, dref

    nfail = 0

    ! (a) G limit of the weight kernel.
    err_max = 0.0_dp
    do k = 1, 50
        tau = real(k, dp) / 51.0_dp
        b = bmin + (bmax - bmin) * tau
        v = 2.5e7_dp
        lam = (2.0_dp * mod(real(k, dp), 2.0_dp) - 1.0_dp) * (0.05_dp + 0.9_dp * tau)
        vpar = lam * v
        mu = 0.5_dp * mass * v**2 * (1.0_dp - lam**2) / b
        hphi = -181.0_dp * (1.0_dp - 0.01_dp * tau)
        rmaj = 160.0_dp + 40.0_dp * tau
        kappa = -2.3e-3_dp * (1.0_dp + tau)
        db1 = 7.0_dp * cos(5.0_dp * tau)
        call evaluate_deltaf_fax_weight(mass, mu, vpar, b, bmax, hphi, rmaj, kappa, &
            db1, -b / (2.0_dp * bmax * lam), dw)
        ref = kappa * mu * (hphi / rmaj)**2 * db1
        err = abs(dw - ref) / abs(ref)
        err_max = max(err_max, err)
    end do
    print '(a, es10.3)', '(a) G limit: max rel error ', err_max
    if (err_max > 1.0e-13_dp) nfail = nfail + 1

    ! (b) coordinates.
    v = 2.5e7_dp
    b = 0.5_dp * (bmin + bmax)
    lam = sqrt(1.0_dp - b / bmax)
    mu = 0.5_dp * mass * v**2 * (1.0_dp - lam**2) / b
    call evaluate_deltaf_fax_coordinates(mass, mu, lam * v, b, bmin, bmax, -1.0_dp, &
        tau, eta, vt)
    print '(a, 3es24.16)', '(b) lambda, eta_hat, vartheta: ', tau, eta, vt
    if (abs(tau - lam) > 1.0e-14_dp .or. abs(eta - 1.0_dp) > 1.0e-14_dp &
        .or. abs(vt + 0.5_dp * pi) > 1.0e-14_dp) nfail = nfail + 1
    call evaluate_deltaf_fax_coordinates(mass, mu, lam * v, bmin, bmin, bmax, 1.0_dp, &
        tau, eta, vt)
    if (abs(vt) > 1.0e-14_dp) nfail = nfail + 1
    call evaluate_deltaf_fax_coordinates(mass, mu, lam * v, bmax, bmin, bmax, 1.0_dp, &
        tau, eta, vt)
    if (abs(vt - pi) > 1.0e-14_dp) nfail = nfail + 1

    ! (c) table: D+ = 1 + 2 vartheta + 3 eta + vartheta eta, D- = -D+.
    call write_table('test_canonical_ntv_fax_table.dat')
    call ntv_fax_load_table('test_canonical_ntv_fax_table.dat')
    err_max = 0.0_dp
    do k = 1, 40
        vt = -pi + 2.0_dp * pi * real(k, dp) / 41.0_dp
        eta = 1.3_dp * real(k, dp) / 41.0_dp
        dref = 1.0_dp + 2.0_dp * vt + 3.0_dp * eta + vt * eta
        err_max = max(err_max, abs(ntv_fax_dchi_deta(vt, eta, .true.) - dref), &
            abs(ntv_fax_dchi_deta(vt, eta, .false.) + dref))
    end do
    d = ntv_fax_dchi_deta(0.0_dp, 5.0_dp, .true.)
    print '(a, es10.3, a, es24.16)', '(c) table max error ', err_max, &
        ', clamped D(0, 5) = ', d
    if (err_max > 1.0e-13_dp .or. abs(d - (1.0_dp + 3.0_dp * 1.3_dp)) > 1.0e-13_dp) &
        nfail = nfail + 1

    if (nfail > 0) then
        print '(a, i0)', 'FAIL: canonical NTV collisional background, failures: ', nfail
        error stop 1
    end if
    print '(a)', 'PASS: canonical NTV collisional background kernels and table'

contains

    subroutine write_table(fname)
        character(len=*), intent(in) :: fname
        integer, parameter :: nt = 9, ne = 5
        real(dp) :: vts(nt), etas(ne), dp_tab(ne, nt)
        integer :: u, i, j

        do j = 1, nt
            vts(j) = -pi + 2.0_dp * pi * real(j - 1, dp) / real(nt - 1, dp)
        end do
        etas = [0.0_dp, 0.4_dp, 0.9_dp, 1.0_dp, 1.3_dp]
        do j = 1, nt
            do i = 1, ne
                dp_tab(i, j) = 1.0_dp + 2.0_dp * vts(j) + 3.0_dp * etas(i) &
                    + vts(j) * etas(i)
            end do
        end do
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(a)') '# test table'
        write(u, '(2i6)') nt, ne
        write(u, '(*(es24.16))') vts
        write(u, '(*(es24.16))') etas
        do j = 1, nt
            write(u, '(*(es24.16))') dp_tab(:, j)
        end do
        do j = 1, nt
            write(u, '(*(es24.16))') -dp_tab(:, j)
        end do
        close(u)
    end subroutine write_table

end program test_canonical_ntv_fax
