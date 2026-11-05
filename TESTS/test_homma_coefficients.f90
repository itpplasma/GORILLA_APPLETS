program test_homma_coefficients

! Unit tests of homma_collision_mod (no orbits):
!   B01  Maxwellian A, D against the oracle of DERIVATIONS/homma_hatayama.wl
!   B02  heat-flux corrections A^(T), D^(T) against the oracle
!   B03  thermal-force moments over test distributions, by quadrature of the
!        pointwise kernel, against the oracle closed forms
!   consistency of A and D: sum_j d(d_ij)/du_j = 2 a_i (from Laplace G = 2 H)
!   continuity at the series/closed-form switch, positivity of D
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

    real(dp), parameter :: tol_oracle = 1.0e-10_dp
    real(dp), parameter :: tol_identity = 1.0e-6_dp
    real(dp), parameter :: tol_switch = 1.0e-11_dp

    character(len=1024) :: ref_dir
    type(homma_units_t) :: units
    logical :: all_passed

    if (command_argument_count() < 1) error stop 'usage: test_homma_coefficients.x <oracle directory>'
    call get_command_argument(1, ref_dir)
    units%k_coulomb = 1.0_dp/(4.0_dp*pi*eps0)
    units%k_lorentz = 1.0_dp

    all_passed = .true.
    call test_coefficient_oracle('B01_oracle.csv', all_passed)
    call test_coefficient_oracle('B02_oracle.csv', all_passed)
    call test_force_moments(all_passed)
    call test_divergence_identity(all_passed)
    call test_series_switch(all_passed)
    call test_positivity(all_passed)

    if (.not. all_passed) then
        print *, 'FAIL: some Homma coefficient tests failed'
        error stop 1
    end if
    print *, 'PASS: all Homma coefficient tests passed'

contains

function hydrogen_background(density, temperature_ev) result(bg)

    real(dp), intent(in) :: density, temperature_ev
    type(homma_background_t) :: bg

    bg%mass = m_p
    bg%charge = e_si
    bg%density = density
    bg%temperature = temperature_ev*e_si
    bg%flow = 0.0_dp
    bg%coulomb_log = 15.0_dp

end function hydrogen_background

subroutine open_oracle(name, unit)

    character(len=*), intent(in) :: name
    integer, intent(out) :: unit
    character(len=4096) :: line

    open(newunit=unit, file=trim(ref_dir)//'/'//name, status='old', action='read')
    do
        read(unit, '(a)') line
        if (line(1:1) /= '#') exit       ! the first non-comment line is the header
    end do

end subroutine open_oracle

function rel_error(computed, reference) result(err)

    ! Largest deviation relative to the largest reference component; exact
    ! agreement is required where the reference vanishes identically.

    real(dp), intent(in) :: computed(:), reference(:)
    real(dp) :: err, scale

    scale = maxval(abs(reference))
    if (scale > 0.0_dp) then
        err = maxval(abs(computed - reference))/scale
    else
        err = merge(0.0_dp, huge(1.0_dp), all(computed == 0.0_dp))
    end if

end function rel_error

function sym6(m) result(s)

    real(dp), intent(in) :: m(3,3)
    real(dp) :: s(6)

    s = [m(1,1), m(1,2), m(1,3), m(2,2), m(2,3), m(3,3)]

end function sym6

subroutine test_coefficient_oracle(name, all_passed)

    ! Rows: Z, gradT(3) [eV/m], x, dir(3), AM(3), AT(3), DM(6), DT(6) in SI.

    character(len=*), intent(in) :: name
    logical, intent(inout) :: all_passed
    type(homma_background_t) :: bg
    type(homma_local_t) :: loc
    real(dp) :: z, grad(3), x, dir(3), ref(18), q(3), v(3)
    real(dp) :: a_m(3), a_t(3), d_m(3,3), d_t(3,3), err(4), max_err(4)
    integer :: unit, ios, nrows

    call open_oracle(name, unit)
    bg = hydrogen_background(1.0e20_dp, 50.0_dp)
    max_err = 0.0_dp
    nrows = 0
    do
        read(unit, *, iostat=ios) z, grad, x, dir, ref
        if (ios /= 0) exit
        nrows = nrows + 1
        q = homma_heat_flux(units, bg, grad*e_si, b_field)
        call homma_local_setup(units, 183.84_dp*amu, z*e_si, bg, q, loc)
        v = x*loc%v_tb*dir
        call homma_coefficients_split(loc, v, a_m, a_t, d_m, d_t)
        err(1) = rel_error(a_m, ref(1:3))
        err(2) = rel_error(a_t, ref(4:6))
        err(3) = rel_error(sym6(d_m), ref(7:12))
        err(4) = rel_error(sym6(d_t), ref(13:18))
        max_err = max(max_err, err)
        if (any(err > tol_oracle)) then
            print *, 'FAIL: ', name, ' row ', nrows, ' x = ', x, ' grad = ', grad, ' errors AM AT DM DT = ', err
            all_passed = .false.
        end if
    end do
    close(unit)

    print '(a,a,i4,a,4es10.2)', name, ': ', nrows, ' rows, max rel. error AM AT DM DT = ', max_err
    if (nrows == 0) then
        print *, 'FAIL: no rows read from ', name
        all_passed = .false.
    end if

end subroutine test_coefficient_oracle

subroutine test_force_moments(all_passed)

    ! F_T = m_a int A^(T) f_a d^3v by quadrature of the pointwise kernel:
    !   maxwellian  Gauss-Hermite in each Cartesian component (T_a = T_b)
    !   shell       Gauss-Legendre in cos(theta) x trapezoid in phi, |u| = x
    !   beam        u = x d
    ! Rows: distribution, x, dir(3), Z, m_a [amu], n [m^-3], T [eV], gradT(3), FT(3).

    logical, intent(inout) :: all_passed
    integer, parameter :: n_gh = 24, n_gl = 16, n_phi = 32
    character(len=32) :: dist
    type(homma_background_t) :: bg
    type(homma_local_t) :: loc
    real(dp) :: x, dir(3), z, ma_amu, dens, t_ev, grad(3), f_ref(3), f(3), q(3)
    real(dp) :: t_gh(n_gh), w_gh(n_gh), t_gl(n_gl), w_gl(n_gl)
    real(dp) :: a_m(3), a_t(3), d_m(3,3), d_t(3,3), v(3), s, ct, st, ph, err, max_err
    integer :: unit, ios, nrows, i, j, k

    call gauss_hermite(n_gh, t_gh, w_gh)
    call gauss_legendre(n_gl, t_gl, w_gl)
    call open_oracle('B03_oracle.csv', unit)
    max_err = 0.0_dp
    nrows = 0
    do
        read(unit, *, iostat=ios) dist, x, dir, z, ma_amu, dens, t_ev, grad, f_ref
        if (ios /= 0) exit
        nrows = nrows + 1
        bg = hydrogen_background(dens, t_ev)
        q = homma_heat_flux(units, bg, grad*e_si, b_field)
        call homma_local_setup(units, ma_amu*amu, z*e_si, bg, q, loc)
        f = 0.0_dp
        select case (trim(dist))
        case ('maxwellian')
            ! v_i = sqrt(2 T_a/m_a) t_i with weight exp(-t^2)/sqrt(pi) per component
            s = sqrt(2.0_dp*t_ev*e_si/(ma_amu*amu))
            do i = 1, n_gh
                do j = 1, n_gh
                    do k = 1, n_gh
                        v = s*[t_gh(i), t_gh(j), t_gh(k)]
                        call homma_coefficients_split(loc, v, a_m, a_t, d_m, d_t)
                        f = f + w_gh(i)*w_gh(j)*w_gh(k)*a_t
                    end do
                end do
            end do
            f = f/pi**1.5_dp
        case ('shell')
            do i = 1, n_gl
                ct = t_gl(i)
                st = sqrt(1.0_dp - ct*ct)
                do j = 1, n_phi
                    ph = 2.0_dp*pi*(j - 1)/n_phi
                    v = x*loc%v_tb*[st*cos(ph), st*sin(ph), ct]
                    call homma_coefficients_split(loc, v, a_m, a_t, d_m, d_t)
                    f = f + w_gl(i)*a_t/(2.0_dp*n_phi)
                end do
            end do
        case ('beam')
            call homma_coefficients_split(loc, x*loc%v_tb*dir, a_m, a_t, d_m, d_t)
            f = a_t
        case default
            print *, 'FAIL: unknown distribution ', trim(dist)
            all_passed = .false.
        end select
        f = ma_amu*amu*f
        err = rel_error(f, f_ref)
        max_err = max(max_err, err)
        if (err > tol_oracle) then
            print *, 'FAIL: B03 row ', nrows, ' ', trim(dist), ' x = ', x, ' Z = ', z, ' error = ', err
            all_passed = .false.
        end if
    end do
    close(unit)

    print '(a,i4,a,es10.2)', 'B03_oracle.csv: ', nrows, ' rows, max rel. error FT = ', max_err
    if (nrows == 0) then
        print *, 'FAIL: no rows read from B03_oracle.csv'
        all_passed = .false.
    end if

end subroutine test_force_moments

subroutine test_divergence_identity(all_passed)

    ! Laplace_u G~ = 2 H~ implies sum_j d(d_ij)/du_j = 2 a_i, separately for
    ! the Maxwellian part and the part linear in g. Central differences. The
    ! error is measured against max(|2 a|, |d|/r): at large r, a^(T) decays
    ! like exp(-r^2) while the terms of d^(T) cancel algebraically.

    logical, intent(inout) :: all_passed
    real(dp), parameter :: h = 1.0e-4_dp
    real(dp) :: u(3), g(3), e(3), div_m(3), div_t(3), err_m, err_t, max_err
    real(dp) :: a_m(3), a_t(3), d_m(3,3), d_t(3,3)
    real(dp) :: a_mp(3), a_tp(3), d_mp(3,3), d_tp(3,3), a_mm(3), a_tm(3), d_mm(3,3), d_tm(3,3)
    real(dp) :: radii(8) = [0.05_dp, 0.3_dp, 0.49_dp, 0.51_dp, 0.8_dp, 1.5_dp, 3.0_dp, 6.0_dp]
    integer :: i, j

    g = [0.3_dp, -0.4_dp, 0.5_dp]
    max_err = 0.0_dp
    do i = 1, size(radii)
        u = radii(i)*[0.48_dp, -0.6_dp, 0.64_dp]
        call homma_kernel_dimensionless(u, g, a_m, a_t, d_m, d_t)
        div_m = 0.0_dp
        div_t = 0.0_dp
        do j = 1, 3
            e = 0.0_dp
            e(j) = h
            call homma_kernel_dimensionless(u + e, g, a_mp, a_tp, d_mp, d_tp)
            call homma_kernel_dimensionless(u - e, g, a_mm, a_tm, d_mm, d_tm)
            div_m = div_m + (d_mp(:,j) - d_mm(:,j))/(2.0_dp*h)
            div_t = div_t + (d_tp(:,j) - d_tm(:,j))/(2.0_dp*h)
        end do
        err_m = maxval(abs(div_m - 2.0_dp*a_m)) &
              /max(2.0_dp*maxval(abs(a_m)), maxval(abs(d_m))/radii(i))
        err_t = maxval(abs(div_t - 2.0_dp*a_t)) &
              /max(2.0_dp*maxval(abs(a_t)), maxval(abs(d_t))/radii(i))
        max_err = max(max_err, err_m, err_t)
        if (max(err_m, err_t) > tol_identity) then
            print *, 'FAIL: div d = 2 a at r = ', radii(i), ' errors (M, T) = ', err_m, err_t
            all_passed = .false.
        end if
    end do
    print '(a,es10.2)', 'identity div d = 2 a: max rel. error = ', max_err

end subroutine test_divergence_identity

subroutine test_series_switch(all_passed)

    logical, intent(inout) :: all_passed
    real(dp) :: below(5), above(5), err

    call homma_kernel_scalars(nearest(homma_r_series, -1.0_dp), below(1), below(2), below(3), below(4), below(5))
    call homma_kernel_scalars(homma_r_series, above(1), above(2), above(3), above(4), above(5))
    err = maxval(abs(above - below)/abs(above))
    print '(a,es10.2)', 'series/closed form at the switch: max rel. jump = ', err
    if (err > tol_switch) then
        print *, 'FAIL: kernel jumps at the series switch: ', below, above
        all_passed = .false.
    end if

end subroutine test_series_switch

subroutine test_positivity(all_passed)

    ! D = D^(M) + D^(T) positive definite for the benchmark gradients over
    ! |u| <= 6 (u along 26 directions); the Monte Carlo relies on it.

    logical, intent(inout) :: all_passed
    real(dp), parameter :: grads(3,4) = reshape([0.0_dp, 0.0_dp, 5.0_dp, 100.0_dp, 0.0_dp, 0.0_dp, &
                                                 300.0_dp, 0.0_dp, 0.0_dp, 300.0_dp, 0.0_dp, 5.0_dp], [3, 4])
    type(homma_background_t) :: bg
    type(homma_local_t) :: loc(1)
    real(dp) :: q(3), a(3), d(3,3), l(3,3), dir(3)
    integer :: ig, ir, i, j, k, ierr, n_fail

    bg = hydrogen_background(1.0e20_dp, 50.0_dp)
    n_fail = 0
    do ig = 1, 4
        q = homma_heat_flux(units, bg, grads(:,ig)*e_si, b_field)
        call homma_local_setup(units, 183.84_dp*amu, 3.0_dp*e_si, bg, q, loc(1))
        do ir = 0, 60
            do i = -1, 1
                do j = -1, 1
                    do k = -1, 1
                        if (i == 0 .and. j == 0 .and. k == 0) cycle
                        dir = real([i, j, k], dp)/norm2(real([i, j, k], dp))
                        call homma_coefficients(loc, 0.1_dp*ir*loc(1)%v_tb*dir, a, d)
                        call homma_cholesky(d, l, ierr)
                        if (ierr /= 0) n_fail = n_fail + 1
                    end do
                end do
            end do
        end do
    end do
    print '(a,i0)', 'positivity of D for benchmark gradients, |u| <= 6: failures = ', n_fail
    if (n_fail > 0) all_passed = .false.

end subroutine test_positivity

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

subroutine gauss_legendre(n, t, w)

    ! Golub-Welsch for weight 1 on [-1, 1].

    integer, intent(in) :: n
    real(dp), intent(out) :: t(n), w(n)
    real(dp) :: offdiag(n-1), zv(n,n), work(2*n)
    integer :: k, info

    t = 0.0_dp
    do k = 1, n - 1
        offdiag(k) = k/sqrt(4.0_dp*k*k - 1.0_dp)
    end do
    call dstev('V', n, t, offdiag, zv, n, work, info)
    if (info /= 0) error stop 'dstev failed in gauss_legendre'
    w = 2.0_dp*zv(1,:)**2

end subroutine gauss_legendre

end program test_homma_coefficients
