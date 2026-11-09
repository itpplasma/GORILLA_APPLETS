program test_homma_kick

! Tests of the stochastic collision step of homma_collision_mod:
!   cholesky       L L^T = 2 D for the benchmark coefficients
!   gc_drift       mean guiding-centre step (noise off), averaged over the
!                  test Maxwellian and the gyrophase by quadrature, against
!                  our first-order B07 targets (DERIVATIONS/homma_hatayama.wl),
!                  for T_b at the guiding centre and on the gyro orbit
!   sampled_kick   sample mean and covariance of dv against A dt, 2 D dt, and
!                  of the guiding-centre shift against the classical-diffusion
!                  covariance R (2 D dt) R^T/Omega^2, R x = -b x x
! W in H+ at n = 1e20 m^-3, T = 50 eV, lnL = 15, B = 1 T e_z, SI units.
! Argument: directory with the oracle files (TESTS/REFERENCE/homma).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use homma_collision_mod

    implicit none

    ! CODATA 2018, SI
    real(dp), parameter :: e_si = 1.602176634e-19_dp
    real(dp), parameter :: eps0 = 8.8541878128e-12_dp
    real(dp), parameter :: amu = 1.66053906660e-27_dp
    real(dp), parameter :: m_p = 1.67262192369e-27_dp
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: b_field(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    real(dp), parameter :: m_w = 183.84_dp*amu

    character(len=1024) :: ref_dir
    type(homma_units_t) :: units
    logical :: all_passed

    if (command_argument_count() < 1) error stop 'usage: test_homma_kick.x <oracle directory>'
    call get_command_argument(1, ref_dir)
    units%k_coulomb = 1.0_dp/(4.0_dp*pi*eps0)
    units%k_lorentz = 1.0_dp

    all_passed = .true.
    call test_cholesky(all_passed)
    call test_gc_drift(all_passed)
    call test_sampled_kick(all_passed)

    if (.not. all_passed) then
        print *, 'FAIL: some Homma kick tests failed'
        error stop 1
    end if
    print *, 'PASS: all Homma kick tests passed'

contains

function benchmark_profile(grad_x_ev_m) result(prof)

    real(dp), intent(in) :: grad_x_ev_m
    type(homma_linear_profile_t) :: prof

    prof%bg0%mass = m_p
    prof%bg0%charge = e_si
    prof%bg0%density = 1.0e20_dp
    prof%bg0%temperature = 50.0_dp*e_si
    prof%bg0%coulomb_log = 15.0_dp
    prof%grad_t = [grad_x_ev_m, 0.0_dp, 0.0_dp]*e_si

end function benchmark_profile

subroutine test_cholesky(all_passed)

    logical, intent(inout) :: all_passed
    type(homma_linear_profile_t) :: prof
    type(homma_local_t) :: loc(1)
    real(dp) :: a(3), d(3,3), l(3,3), v(3), err, max_err, bad(3,3)
    integer :: i, ierr

    prof = benchmark_profile(300.0_dp)
    prof%grad_t(3) = 5.0_dp*e_si
    call homma_local_at(units, m_w, 3.0_dp*e_si, prof, [0.0_dp, 0.0_dp, 0.0_dp], b_field, loc(1), ierr)
    max_err = 0.0_dp
    do i = 1, 50
        v = loc(1)%v_tb*0.1_dp*i*[sin(0.7_dp*i)*cos(1.3_dp*i), sin(0.7_dp*i)*sin(1.3_dp*i), cos(0.7_dp*i)]
        call homma_coefficients(loc, v, a, d)
        call homma_cholesky(2.0_dp*d, l, ierr)
        if (ierr /= 0) then
            print *, 'FAIL: Cholesky failed at |u| = ', 0.1_dp*i
            all_passed = .false.
            cycle
        end if
        err = maxval(abs(matmul(l, transpose(l)) - 2.0_dp*d))/maxval(abs(2.0_dp*d))
        max_err = max(max_err, err)
    end do
    print '(a,es10.2)', 'cholesky: max rel. error of L L^T - 2D = ', max_err
    if (max_err > 1.0e-14_dp) all_passed = .false.

    ! a matrix that is not positive definite must be reported
    bad = reshape([1.0_dp, 2.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
    call homma_cholesky(bad, l, ierr)
    if (ierr /= 1) then
        print *, 'FAIL: indefinite matrix not reported'
        all_passed = .false.
    end if

end subroutine test_cholesky

subroutine test_gc_drift(all_passed)

    ! Mean step with xi = 0, i.e. the expectation over the Wiener increment
    ! (dx_gc is linear in xi). The test Maxwellian at T_a = 50 eV is
    ! integrated by Gauss-Hermite in (v_x, v_y, v_z); (v_perp, phase) follow
    ! from the basis used by homma_gc_step. dt is arbitrary for the mean.
    ! Targets are first order in grad T. The guiding-centre mode is linear in
    ! grad T and must match to quadrature precision. The gyro-orbit mode is
    ! not: the kernel at T(X + rho) adds a relative correction of second
    ! order, about 1.2e-3 (grad T/100 eV/m)^2 for W1+ (rho = 10 mm), less for
    ! higher charge states. The tolerance scales with (grad T)^2, so 1 eV/m
    ! checks the first-order part to 2e-7.

    logical, intent(inout) :: all_passed
    integer, parameter :: n_gh = 20
    character(len=4096) :: line
    real(dp) :: t_gh(n_gh), w_gh(n_gh), e1(3), e2(3), target(6), z, grad_x
    real(dp) :: mean_gc(3), mean_orbit(3), scale_grad
    integer :: unit, ios, igrad
    real(dp), parameter :: tol(2) = [2.0e-3_dp, 2.0e-7_dp]   ! 100 eV/m, 1 eV/m (gyro orbit)

    call gauss_hermite(n_gh, t_gh, w_gh)
    call homma_perp_basis(b_field, e1, e2)

    open(newunit=unit, file=trim(ref_dir)//'/B07_targets.csv', status='old', action='read')
    do
        read(unit, '(a)') line
        if (line(1:1) /= '#') exit
    end do
    do
        read(unit, *, iostat=ios) z, grad_x, target
        if (ios /= 0) exit
        if (grad_x /= 100.0_dp) cycle
        do igrad = 1, 2
            scale_grad = merge(1.0_dp, 0.01_dp, igrad == 1)
            call mean_drift(t_gh, w_gh, e1, e2, z, scale_grad*grad_x, homma_at_guiding_centre, mean_gc)
            call mean_drift(t_gh, w_gh, e1, e2, z, scale_grad*grad_x, homma_at_gyro_orbit, mean_orbit)
            print '(a,i3,a,f6.1,a,2f11.6,a,f11.6,a,2f11.6,a,f11.6)', 'gc_drift Z =', nint(z), &
                ', grad T =', scale_grad*grad_x, ' eV/m: guiding centre dX, dY =', mean_gc(1:2), &
                ' (target', scale_grad*target(1), '), gyro orbit dX, dY =', mean_orbit(1:2), &
                ' (target', scale_grad*target(4), ')'
            ! guiding-centre mode is linear in grad T: quadrature precision
            call compare('guiding centre dX/dt', mean_gc(1), scale_grad*target(1), 1.0e-9_dp, all_passed)
            call compare('guiding centre dY/dt', mean_gc(2), scale_grad*target(2), 1.0e-6_dp, all_passed)
            call compare('gyro orbit dX/dt', mean_orbit(1), scale_grad*target(4), tol(igrad), all_passed)
            if (abs(mean_gc(3)) + abs(mean_orbit(3)) > 1.0e-12_dp*abs(target(1))) then
                print *, 'FAIL: parallel guiding-centre shift ', mean_gc(3), mean_orbit(3)
                all_passed = .false.
            end if
        end do
    end do
    close(unit)


end subroutine test_gc_drift

subroutine mean_drift(t_gh, w_gh, e1, e2, z, grad_x_ev_m, background_at, mean)

    ! Maxwellian (T_a = 50 eV) and gyrophase average of the mean step dx_gc/dt.

    real(dp), intent(in) :: t_gh(:), w_gh(:), e1(3), e2(3), z, grad_x_ev_m
    real(dp), parameter :: dt = 1.0e-7_dp
    integer, intent(in) :: background_at
    real(dp), intent(out) :: mean(3)
    type(homma_linear_profile_t) :: prof
    real(dp) :: s, v(3), vpar, vperp, phase_u, dx(3), w
    integer :: i, j, k, ierr

    prof = benchmark_profile(grad_x_ev_m)
    s = sqrt(2.0_dp*50.0_dp*e_si/m_w)
    mean = 0.0_dp
    do i = 1, size(t_gh)
        do j = 1, size(t_gh)
            do k = 1, size(t_gh)
                v = s*[t_gh(i), t_gh(j), t_gh(k)]
                w = w_gh(i)*w_gh(j)*w_gh(k)/pi**1.5_dp
                vpar = dot_product(v, b_field)
                vperp = norm2(v - vpar*b_field)
                phase_u = modulo(atan2(dot_product(v, e2), dot_product(v, e1))/(2.0_dp*pi), 1.0_dp)
                call homma_gc_step(units, m_w, z*e_si, prof, b_field, background_at, &
                                   [0.0_dp, 0.0_dp, 0.0_dp], vpar, vperp, dt, phase_u, &
                                   [0.0_dp, 0.0_dp, 0.0_dp], dx, ierr)
                if (ierr /= 0) error stop 'homma_gc_step failed in test_gc_drift'
                mean = mean + w*dx/dt
            end do
        end do
    end do

end subroutine mean_drift


subroutine test_sampled_kick(all_passed)

    ! Fixed velocity, random gyrophase and Wiener increments. For the
    ! velocity kick the gyrophase is fixed by v; for the guiding-centre shift
    ! the covariance is checked at that phase as well.

    logical, intent(inout) :: all_passed
    integer, parameter :: n = 400000
    real(dp), parameter :: dt = 1.0e-7_dp
    type(homma_linear_profile_t) :: prof
    type(homma_local_t) :: loc(1)
    real(dp) :: v(3), a(3), d(3,3), dv(3), xi(3), sum1(3), sum2(3,3), mean(3), cov(3,3)
    real(dp) :: sx1(3), sx2(3,3), covx(3,3), r(3,3), expect_x(3,3), omega, se(3)
    real(dp) :: vpar, vperp, vpar0, vperp0, phase_u, dx(3), e1(3), e2(3), err
    integer :: i, ierr, seed_size
    integer, allocatable :: seed(:)

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 20261109
    call random_seed(put=seed)

    prof = benchmark_profile(300.0_dp)
    prof%grad_t(3) = 5.0_dp*e_si
    call homma_local_at(units, m_w, 3.0_dp*e_si, prof, [0.0_dp, 0.0_dp, 0.0_dp], b_field, loc(1), ierr)

    ! v at |u| = 0.5, pitch 0.6, gyrophase 1 rad in the step's basis
    call homma_perp_basis(b_field, e1, e2)
    vpar0 = 0.5_dp*loc(1)%v_tb*0.6_dp
    vperp0 = 0.5_dp*loc(1)%v_tb*0.8_dp
    phase_u = 1.0_dp/(2.0_dp*pi)
    v = vpar0*b_field + vperp0*(cos(1.0_dp)*e1 + sin(1.0_dp)*e2)
    call homma_coefficients(loc, v, a, d)

    sum1 = 0.0_dp; sum2 = 0.0_dp; sx1 = 0.0_dp; sx2 = 0.0_dp
    do i = 1, n
        call normal3(xi)
        call homma_kick(loc, v, dt, xi, dv, ierr)
        sum1 = sum1 + dv
        sum2 = sum2 + spread(dv, 2, 3)*spread(dv, 1, 3)
        vpar = vpar0
        vperp = vperp0
        call homma_gc_step(units, m_w, 3.0_dp*e_si, prof, b_field, homma_at_guiding_centre, &
                           [0.0_dp, 0.0_dp, 0.0_dp], vpar, vperp, dt, phase_u, xi, dx, ierr)
        sx1 = sx1 + dx
        sx2 = sx2 + spread(dx, 2, 3)*spread(dx, 1, 3)
    end do
    mean = sum1/n
    cov = sum2/n - spread(mean, 2, 3)*spread(mean, 1, 3)
    covx = sx2/n - spread(sx1/n, 2, 3)*spread(sx1/n, 1, 3)

    ! mean within 5 standard errors, covariances to 1 % (statistical 0.2 %)
    se = sqrt(2.0_dp*[d(1,1), d(2,2), d(3,3)]*dt/n)
    print '(a,3f8.2)', 'sampled_kick: (mean dv - A dt)/standard error = ', (mean - a*dt)/se
    if (any(abs(mean - a*dt) > 5.0_dp*se)) all_passed = .false.

    err = maxval(abs(cov - 2.0_dp*d*dt))/maxval(abs(2.0_dp*d*dt))
    print '(a,es10.2)', 'sampled_kick: rel. error of cov(dv) - 2 D dt = ', err
    if (err > 1.0e-2_dp) all_passed = .false.

    omega = 3.0_dp*e_si*norm2(b_field)/m_w
    r = reshape([0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [3, 3])  ! x -> -e_z x x
    expect_x = matmul(r, matmul(2.0_dp*d*dt, transpose(r)))/omega**2
    err = maxval(abs(covx - expect_x))/maxval(abs(expect_x))
    print '(a,es10.2)', 'sampled_kick: rel. error of cov(dX) - R 2D dt R^T/Omega^2 = ', err
    if (err > 1.0e-2_dp) all_passed = .false.

end subroutine test_sampled_kick

subroutine normal3(xi)

    real(dp), intent(out) :: xi(3)
    real(dp) :: u(4)

    call random_number(u)
    u(1) = max(u(1), tiny(1.0_dp))
    u(3) = max(u(3), tiny(1.0_dp))
    xi(1) = sqrt(-2.0_dp*log(u(1)))*cos(2.0_dp*pi*u(2))
    xi(2) = sqrt(-2.0_dp*log(u(1)))*sin(2.0_dp*pi*u(2))
    xi(3) = sqrt(-2.0_dp*log(u(3)))*cos(2.0_dp*pi*u(4))

end subroutine normal3

subroutine compare(label, value, reference, tol, all_passed)

    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, reference, tol
    logical, intent(inout) :: all_passed

    if (abs(value - reference) > tol*abs(reference)) then
        print '(a,a,a,es22.14,a,es22.14,a,es10.2)', 'FAIL: ', label, ' = ', value, ', target ', reference, &
            ', rel. error ', abs(value - reference)/abs(reference)
        all_passed = .false.
    end if

end subroutine compare

subroutine gauss_hermite(n, t, w)

    ! Golub-Welsch: nodes and weights for weight exp(-t^2) on the real line.

    integer, intent(in) :: n
    real(dp), intent(out) :: t(n), w(n)
    real(dp) :: offdiag(n-1), zv(n,n), work(2*n)
    integer :: k, info

    t = 0.0_dp
    do k = 1, n - 1
        offdiag(k) = sqrt(0.5_dp*k)
    end do
    call dstev('V', n, t, offdiag, zv, n, work, info)
    if (info /= 0) error stop 'dstev failed in gauss_hermite'
    w = sqrt(pi)*zv(1,:)**2

end subroutine gauss_hermite

end program test_homma_kick
