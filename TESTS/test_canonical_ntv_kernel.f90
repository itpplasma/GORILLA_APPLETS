!> Gate (b) of the canonical NTV estimator: the generated torque sample
!> -dH1/dphi equals the finite-difference phi derivative of the guiding-center
!> Hamiltonian H(phi) = m vpar^2/2 + mu |B| + e Phi1 at fixed canonical
!> momentum p_phi, magnetic moment mu, flux label and poloidal angle, sampled
!> along a test orbit through the analytic circular field with one harmonic.
!> vpar follows from p_phi = e A_phi/c + m vpar B_phi/|B| with axisymmetric A_phi
!> and covariant B_phi, so the field-direction term m vpar^2/B dB1/dphi enters.
!> The harmonic goes through canonical_ntv_mod (the estimator's own path).
program test_canonical_ntv_kernel

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use deltaf_circular_chart_kernel, only: evaluate_deltaf_circular_chart
    use deltaf_weight_kernel, only: evaluate_deltaf_weight
    use canonical_ntv_mod, only: ntv_m_pol, ntv_n_tor, ntv_amplitude, ntv_bmn, &
        ntv_phimn, ntv_harmonic_dphi

    implicit none

    real(dp), parameter :: R0 = 165.0_dp, a = 70.0_dp, B0 = 18000.0_dp
    real(dp), parameter :: q0 = 1.0_dp, q1 = 2.0_dp
    real(dp), parameter :: mass = 3.3436e-24_dp, charge = 4.8032e-10_dp
    real(dp), parameter :: c_light = 2.9979e10_dp, pi = acos(-1.0_dp)
    real(dp), parameter :: dphi = 1.0e-4_dp
    integer, parameter :: n_samples = 64
    real(dp) :: tau, r, theta, phi, vpar0, mu, pphi, bphi_cov, aphi
    real(dp) :: bmod0, s, h_phi, dum(8), db1_dphi, dp1_dphi, B, vpar
    real(dp) :: ln_f0, dlnf0, dh1, dwdt, tq, fd, hp, hm, err, err_max, scale
    real(dp) :: err_no_mirror, err_no_mirror_min
    integer :: k

    ntv_m_pol = 3
    ntv_n_tor = 2
    ntv_amplitude = 1.0_dp
    ntv_bmn = 1.0e-2_dp * B0
    ntv_phimn = 50.0_dp / 299.79_dp

    err_max = 0.0_dp
    err_no_mirror_min = huge(1.0_dp)
    do k = 1, n_samples
        tau = 2.0_dp * pi * real(k - 1, dp) / real(n_samples, dp)
        r = a * (0.5_dp + 0.2_dp * sin(tau))
        theta = 3.0_dp * tau
        phi = 0.7_dp * tau
        vpar0 = 1.2e8_dp * cos(2.0_dp * tau)
        call evaluate_deltaf_circular_chart(R0, a, B0, q0, q1, r, theta, bmod0, &
            dum(1), dum(2), aphi, dum(3), dum(4), dum(5), h_phi, s, dum(6), dum(7))
        bphi_cov = h_phi * bmod0
        mu = 0.5_dp * mass * (0.9e8_dp)**2 / bmod0
        call ntv_harmonic_dphi(s, theta, phi, db1_dphi, dp1_dphi)
        B = bmod0 + field_pert(s, theta, phi)
        vpar = vpar0
        pphi = charge * aphi / c_light + mass * vpar * bphi_cov / B

        hp = hamiltonian(phi + dphi)
        hm = hamiltonian(phi - dphi)
        fd = (hp - hm) / (2.0_dp * dphi)

        call evaluate_deltaf_weight(mass, charge, c_light, mu, vpar, B, 1.0e-9_dp, &
            1.0e13_dp, -1.0e5_dp, 1.6e-9_dp, -1.0e-17_dp, 0.0_dp, 0.0_dp, db1_dphi, &
            dp1_dphi, ln_f0, dlnf0, dh1, dwdt, tq)
        scale = (mu + mass * vpar**2 / B) * abs(db1_dphi) + charge * abs(dp1_dphi)
        err = abs(-tq - fd) / scale
        err_max = max(err_max, err)
        ! Without the field-direction term the comparison must fail clearly.
        err_no_mirror = abs(mu * db1_dphi + charge * dp1_dphi - fd) / scale
        if (abs(vpar) > 0.3e8_dp) err_no_mirror_min = min(err_no_mirror_min, err_no_mirror)
        if (abs(tq + dh1) > 0.0_dp) error stop 'torque sample is not -dH1/dphi'
        if (abs(dwdt - dlnf0 * dh1) > 1.0e-14_dp * abs(dwdt) + tiny(1.0_dp)) &
            error stop 'dw/dt is not dlnf0/dpphi dH1/dphi'
    end do
    print '(a, es10.3, a, es10.3)', 'max relative |(-torque) - FD dH/dphi| = ', &
        err_max, '; without the m vpar^2/B term (min) = ', err_no_mirror_min
    if (err_max > 1.0e-6_dp) error stop 'generated dH1/dphi disagrees with FD'
    if (err_no_mirror_min < 1.0e-2_dp) error stop 'FD oracle insensitive to mirror term'
    print '(a)', 'PASS: generated torque sample equals FD dH/dphi along the orbit'

contains

    real(dp) function field_pert(s_l, th, ph) result(db1)
        real(dp), intent(in) :: s_l, th, ph
        db1 = ntv_amplitude * ntv_bmn * cos(real(ntv_m_pol, dp) * th &
                                            - real(ntv_n_tor, dp) * ph)
    end function field_pert

    real(dp) function hamiltonian(ph) result(H)
        real(dp), intent(in) :: ph
        real(dp) :: Bl, vp
        Bl = bmod0 + field_pert(s, theta, ph)
        vp = (pphi - charge * aphi / c_light) * Bl / (mass * bphi_cov)
        H = 0.5_dp * mass * vp**2 + mu * Bl + charge * ntv_amplitude * ntv_phimn &
            * cos(real(ntv_m_pol, dp) * theta - real(ntv_n_tor, dp) * ph)
    end function hamiltonian

end program test_canonical_ntv_kernel
