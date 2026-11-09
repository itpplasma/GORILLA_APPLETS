module homma_collision_mod

! Fokker-Planck collision coefficients of a test ion in background ions whose
! Maxwellian is distorted by the classical heat flux (Homma & Hatayama 2014,
! Contrib. Plasma Phys. 54, 394, Eqs. 1, 2, 5, 6), and the Euler-Maruyama
! velocity kick dv = A dt + L xi sqrt(dt) with L L^T = 2 D (Ito).
!
! Derivation and oracle values: DERIVATIONS/homma_hatayama.wl
! Design: docs/plans/2026-11-05-homma-thermal-force-design.md
!
! The module is independent of the GORILLA state and of the unit system. A
! homma_units_t holds the two factors in which SI and Gaussian units differ:
!   k_coulomb = 1/(4 pi eps0) (SI) or 1 (Gaussian)   in q_a q_b k_coulomb/r
!   k_lorentz = 1 (SI) or 1/c (Gaussian)              in Omega = k_lorentz q B/m
! Temperatures are energies (J or erg).
!
! Normalisation per background species b:
!   v_tb = sqrt(2 T_b/m_b),  u = (v - u_b)/v_tb,  r = |u|,
!   g = m_b v_tb q/(n_b T_b^2)                    (dimensionless heat flux)
!   A = p_a a(u, g),   p_a = Gamma (1 + m_a/m_b) n_b/v_tb^2
!   D = p_d d(u, g),   p_d = Gamma n_b/(2 v_tb)
!   Gamma = 4 pi k_coulomb^2 q_a^2 q_b^2 lnL/m_a^2
!   a = -alpha(r) u - k_t exp(-r^2) (g - 2 (g.u) u),        k_t = 2/(5 sqrt(pi))
!   d = d1(r) I + d2(r) u u + t1(r) (g u + u g + (g.u) I) + t2(r) (g.u) u u

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    include 'homma_series_coefficients.inc'

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    real(dp), parameter :: sqrt_pi = 1.77245385090551602729816748334114518_dp
    real(dp), parameter :: k_thermal = 2.0_dp/(5.0_dp*sqrt_pi)

    type, public :: homma_units_t
        real(dp) :: k_coulomb = 1.0_dp
        real(dp) :: k_lorentz = 1.0_dp
    end type homma_units_t

    type, public :: homma_background_t
        real(dp) :: mass = 0.0_dp
        real(dp) :: charge = 0.0_dp
        real(dp) :: density = 0.0_dp
        real(dp) :: temperature = 0.0_dp       ! energy units
        real(dp) :: flow(3) = 0.0_dp
        real(dp) :: coulomb_log = 0.0_dp
    end type homma_background_t

    ! Prefactors of one test species colliding with one background species at
    ! one position, from homma_local_setup.
    type, public :: homma_local_t
        real(dp) :: v_tb = 0.0_dp
        real(dp) :: p_a = 0.0_dp
        real(dp) :: p_d = 0.0_dp
        real(dp) :: g(3) = 0.0_dp
        real(dp) :: flow(3) = 0.0_dp
    end type homma_local_t

    ! Background with a constant temperature gradient (benchmark profiles):
    ! T(x) = bg0%temperature + grad_t.(x - x0), other quantities uniform. The
    ! heat flux is evaluated with the local T(x).
    type, public :: homma_linear_profile_t
        type(homma_background_t) :: bg0
        real(dp) :: x0(3) = 0.0_dp
        real(dp) :: grad_t(3) = 0.0_dp          ! energy per length
        logical :: heat_flux = .true.           ! .false.: Maxwellian operator only
    end type homma_linear_profile_t

    ! Where homma_gc_step evaluates the background
    integer, parameter, public :: homma_at_guiding_centre = 1
    integer, parameter, public :: homma_at_gyro_orbit = 2

    public :: homma_heat_flux, homma_local_setup, homma_kernel_scalars, homma_kernel_dimensionless
    public :: homma_coefficients_split, homma_coefficients, homma_cholesky, homma_kick
    public :: homma_local_at, homma_perp_basis, homma_gc_step
    public :: homma_n_series, homma_r_series

contains

pure function homma_heat_flux(units, bg, grad_t, b_field) result(q)

    ! Classical ion heat flux printed below Homma Eq. (1) (Braginskii):
    !   q = -kappa_par grad_par T + kappa_wedge b x grad_perp T - kappa_perp grad_perp T
    ! with tau_b as printed by Homma, 12 pi^(3/2) eps0^2 sqrt(m_b) T^(3/2)/(n q_b^4 lnL).
    ! grad_t in energy per length, b_field the magnetic field vector.

    type(homma_units_t), intent(in) :: units
    type(homma_background_t), intent(in) :: bg
    real(dp), intent(in) :: grad_t(3), b_field(3)
    real(dp) :: q(3)
    real(dp) :: b_abs, b_hat(3), omega_b, tau_b, kappa_par, kappa_wedge, kappa_perp
    real(dp) :: grad_par(3), grad_perp(3)

    b_abs = norm2(b_field)
    b_hat = b_field/b_abs
    omega_b = units%k_lorentz*abs(bg%charge)*b_abs/bg%mass
    tau_b = 3.0_dp/(4.0_dp*sqrt_pi)*sqrt(bg%mass)*bg%temperature**1.5_dp &
          /(units%k_coulomb**2*bg%density*bg%charge**4*bg%coulomb_log)

    kappa_par = 3.9_dp*bg%density*bg%temperature*tau_b/bg%mass
    kappa_wedge = 5.0_dp*bg%density*bg%temperature/(2.0_dp*bg%mass*omega_b)
    kappa_perp = bg%density*bg%temperature/(bg%mass*omega_b**2*tau_b)

    grad_par = dot_product(grad_t, b_hat)*b_hat
    grad_perp = grad_t - grad_par
    q = -kappa_par*grad_par + kappa_wedge*cross(b_hat, grad_perp) - kappa_perp*grad_perp

end function homma_heat_flux

pure subroutine homma_local_setup(units, test_mass, test_charge, bg, q, loc)

    ! Prefactors for test species (test_mass, test_charge) in background bg
    ! with heat flux q. q = 0 gives the Maxwellian operator.

    type(homma_units_t), intent(in) :: units
    real(dp), intent(in) :: test_mass, test_charge
    type(homma_background_t), intent(in) :: bg
    real(dp), intent(in) :: q(3)
    type(homma_local_t), intent(out) :: loc
    real(dp) :: gamma_ab

    gamma_ab = 4.0_dp*pi*units%k_coulomb**2*test_charge**2*bg%charge**2*bg%coulomb_log/test_mass**2
    loc%v_tb = sqrt(2.0_dp*bg%temperature/bg%mass)
    loc%p_a = gamma_ab*(1.0_dp + test_mass/bg%mass)*bg%density/loc%v_tb**2
    loc%p_d = gamma_ab*bg%density/(2.0_dp*loc%v_tb)
    loc%g = bg%mass*loc%v_tb*q/(bg%density*bg%temperature**2)
    loc%flow = bg%flow

end subroutine homma_local_setup

pure subroutine homma_kernel_scalars(r, alpha, d1, d2, t1, t2)

    ! Scalar kernel functions of r = |u|: Taylor series in w = r^2 below
    ! homma_r_series, where the closed forms cancel, closed forms above.

    real(dp), intent(in) :: r
    real(dp), intent(out) :: alpha, d1, d2, t1, t2
    real(dp) :: w, er, phi_p, r2, r3, r4, r5, r6, r7

    if (r < homma_r_series) then
        w = r*r
        alpha = horner(c_alpha, w)
        d1 = horner(c_d1, w)
        d2 = horner(c_d2, w)
        t1 = horner(c_t1, w)
        t2 = horner(c_t2, w)
        return
    end if

    er = erf(r)
    phi_p = 2.0_dp/sqrt_pi*exp(-r*r)
    r2 = r*r; r3 = r2*r; r4 = r2*r2; r5 = r4*r; r6 = r3*r3; r7 = r6*r

    alpha = (er - r*phi_p)/r3
    d1 = 0.5_dp*((2.0_dp/r - 1.0_dp/r3)*er + phi_p/r2)
    d2 = 0.5_dp*((-2.0_dp/r3 + 3.0_dp/r5)*er - 3.0_dp/r4*phi_p)
    t1 = 0.2_dp*(-1.5_dp/r5*er + (1.0_dp/r2 + 1.5_dp/r4)*phi_p)
    t2 = 0.2_dp*(7.5_dp/r7*er + (-2.0_dp/r2 - 5.0_dp/r4 - 7.5_dp/r6)*phi_p)

end subroutine homma_kernel_scalars

pure subroutine homma_kernel_dimensionless(u, g, a_m, a_t, d_m, d_t)

    ! Dimensionless a = grad_u H~ and d = grad_u grad_u G~, split into the
    ! Maxwellian part (m) and the part linear in g (t).

    real(dp), intent(in) :: u(3), g(3)
    real(dp), intent(out) :: a_m(3), a_t(3), d_m(3,3), d_t(3,3)
    real(dp) :: r, alpha, d1, d2, t1, t2, gu, uu(3,3)
    integer :: i

    r = norm2(u)
    call homma_kernel_scalars(r, alpha, d1, d2, t1, t2)
    gu = dot_product(g, u)
    uu = spread(u, 2, 3)*spread(u, 1, 3)

    a_m = -alpha*u
    a_t = -k_thermal*exp(-r*r)*(g - 2.0_dp*gu*u)

    d_m = d2*uu
    d_t = t1*(spread(g, 2, 3)*spread(u, 1, 3) + spread(u, 2, 3)*spread(g, 1, 3)) + t2*gu*uu
    do i = 1, 3
        d_m(i,i) = d_m(i,i) + d1
        d_t(i,i) = d_t(i,i) + t1*gu
    end do

end subroutine homma_kernel_dimensionless

pure subroutine homma_coefficients_split(loc, v, a_m, a_t, d_m, d_t)

    ! Physical A and D of one background species at test velocity v,
    ! Maxwellian part (m) and heat-flux part (t).

    type(homma_local_t), intent(in) :: loc
    real(dp), intent(in) :: v(3)
    real(dp), intent(out) :: a_m(3), a_t(3), d_m(3,3), d_t(3,3)

    call homma_kernel_dimensionless((v - loc%flow)/loc%v_tb, loc%g, a_m, a_t, d_m, d_t)
    a_m = loc%p_a*a_m
    a_t = loc%p_a*a_t
    d_m = loc%p_d*d_m
    d_t = loc%p_d*d_t

end subroutine homma_coefficients_split

pure subroutine homma_coefficients(locs, v, a, d)

    ! Total A and D summed over the background species in locs.

    type(homma_local_t), intent(in) :: locs(:)
    real(dp), intent(in) :: v(3)
    real(dp), intent(out) :: a(3), d(3,3)
    real(dp) :: a_m(3), a_t(3), d_m(3,3), d_t(3,3)
    integer :: k

    a = 0.0_dp
    d = 0.0_dp
    do k = 1, size(locs)
        call homma_coefficients_split(locs(k), v, a_m, a_t, d_m, d_t)
        a = a + a_m + a_t
        d = d + d_m + d_t
    end do

end subroutine homma_coefficients

pure subroutine homma_cholesky(m, l, ierr)

    ! Lower-triangular l with l l^T = m for a symmetric 3x3 m.
    ! ierr = 1 if m is not positive definite; l is then undefined.

    real(dp), intent(in) :: m(3,3)
    real(dp), intent(out) :: l(3,3)
    integer, intent(out) :: ierr
    real(dp) :: piv

    ierr = 0
    l = 0.0_dp

    piv = m(1,1)
    if (piv <= 0.0_dp) then
        ierr = 1
        return
    end if
    l(1,1) = sqrt(piv)
    l(2,1) = m(2,1)/l(1,1)
    l(3,1) = m(3,1)/l(1,1)

    piv = m(2,2) - l(2,1)**2
    if (piv <= 0.0_dp) then
        ierr = 1
        return
    end if
    l(2,2) = sqrt(piv)
    l(3,2) = (m(3,2) - l(3,1)*l(2,1))/l(2,2)

    piv = m(3,3) - l(3,1)**2 - l(3,2)**2
    if (piv <= 0.0_dp) then
        ierr = 1
        return
    end if
    l(3,3) = sqrt(piv)

end subroutine homma_cholesky

pure subroutine homma_kick(locs, v, dt, xi, dv, ierr)

    ! Euler-Maruyama step of dv = A dt + L dW, L L^T = 2 D, <dW dW^T> = I dt.
    ! xi: three independent standard normal numbers.
    ! ierr = 1 if D is not positive definite; dv is then zero and nothing is
    ! repaired, the caller counts and reports it.

    type(homma_local_t), intent(in) :: locs(:)
    real(dp), intent(in) :: v(3), dt, xi(3)
    real(dp), intent(out) :: dv(3)
    integer, intent(out) :: ierr
    real(dp) :: a(3), d(3,3), l(3,3)

    call homma_coefficients(locs, v, a, d)
    call homma_cholesky(2.0_dp*d, l, ierr)
    if (ierr /= 0) then
        dv = 0.0_dp
        return
    end if
    dv = a*dt + matmul(l, xi)*sqrt(dt)

end subroutine homma_kick

pure subroutine homma_local_at(units, test_mass, test_charge, prof, x, b_field, loc, ierr)

    ! Prefactors at position x for a linear temperature profile.
    ! ierr = 2 if T(x) <= 0.

    type(homma_units_t), intent(in) :: units
    real(dp), intent(in) :: test_mass, test_charge
    type(homma_linear_profile_t), intent(in) :: prof
    real(dp), intent(in) :: x(3), b_field(3)
    type(homma_local_t), intent(out) :: loc
    integer, intent(out) :: ierr
    type(homma_background_t) :: bg
    real(dp) :: q(3)

    ierr = 0
    bg = prof%bg0
    bg%temperature = prof%bg0%temperature + dot_product(prof%grad_t, x - prof%x0)
    if (bg%temperature <= 0.0_dp) then
        ierr = 2
        return
    end if
    q = 0.0_dp
    if (prof%heat_flux) q = homma_heat_flux(units, bg, prof%grad_t, b_field)
    call homma_local_setup(units, test_mass, test_charge, bg, q, loc)

end subroutine homma_local_at

pure subroutine homma_perp_basis(b_hat, e1, e2)

    ! Unit vectors e1, e2 with (e1, e2, b_hat) right-handed and orthonormal.

    real(dp), intent(in) :: b_hat(3)
    real(dp), intent(out) :: e1(3), e2(3)
    real(dp) :: a(3)

    if (abs(b_hat(1)) < 0.9_dp) then
        a = [1.0_dp, 0.0_dp, 0.0_dp]
    else
        a = [0.0_dp, 1.0_dp, 0.0_dp]
    end if
    e1 = a - dot_product(a, b_hat)*b_hat
    e1 = e1/norm2(e1)
    e2 = cross(b_hat, e1)

end subroutine homma_perp_basis

pure subroutine homma_gc_step(units, test_mass, test_charge, prof, b_field, background_at, &
                              x_gc, vpar, vperp, dt, phase_u, xi, dx_gc, ierr)

    ! One collision step of a guiding centre at x_gc with (vpar, vperp):
    !   1. gyrophase phi = 2 pi phase_u, v = vpar b + vperp (cos phi e1 + sin phi e2)
    !   2. background at x_gc (homma_at_guiding_centre) or at the particle
    !      x_gc + rho, rho = b x v/Omega_a (homma_at_gyro_orbit)
    !   3. Euler-Maruyama kick dv (homma_kick)
    !   4. vpar, vperp from v + dv; guiding-centre shift dx_gc = -b x dv/Omega_a,
    !      the particle position stays fixed.
    ! Omega_a = k_lorentz q_a |B|/m_a carries the sign of the charge.
    ! phase_u: uniform on [0, 1); xi: three standard normal numbers.
    ! ierr /= 0 (from homma_local_at or homma_kick): vpar, vperp unchanged, dx_gc = 0.

    type(homma_units_t), intent(in) :: units
    real(dp), intent(in) :: test_mass, test_charge
    type(homma_linear_profile_t), intent(in) :: prof
    real(dp), intent(in) :: b_field(3), x_gc(3), dt, phase_u, xi(3)
    integer, intent(in) :: background_at
    real(dp), intent(inout) :: vpar, vperp
    real(dp), intent(out) :: dx_gc(3)
    integer, intent(out) :: ierr
    type(homma_local_t) :: loc(1)
    real(dp) :: b_abs, b_hat(3), e1(3), e2(3), omega, phi, v(3), dv(3), x_eval(3)

    dx_gc = 0.0_dp
    b_abs = norm2(b_field)
    b_hat = b_field/b_abs
    omega = units%k_lorentz*test_charge*b_abs/test_mass

    call homma_perp_basis(b_hat, e1, e2)
    phi = 2.0_dp*pi*phase_u
    v = vpar*b_hat + vperp*(cos(phi)*e1 + sin(phi)*e2)

    x_eval = x_gc
    if (background_at == homma_at_gyro_orbit) x_eval = x_gc + cross(b_hat, v)/omega

    call homma_local_at(units, test_mass, test_charge, prof, x_eval, b_field, loc(1), ierr)
    if (ierr /= 0) return
    call homma_kick(loc, v, dt, xi, dv, ierr)
    if (ierr /= 0) return

    v = v + dv
    vpar = dot_product(v, b_hat)
    vperp = norm2(v - vpar*b_hat)
    dx_gc = -cross(b_hat, dv)/omega

end subroutine homma_gc_step

pure function horner(c, w) result(s)

    real(dp), intent(in) :: c(0:), w
    real(dp) :: s
    integer :: k

    s = c(ubound(c, 1))
    do k = ubound(c, 1) - 1, 0, -1
        s = s*w + c(k)
    end do

end function horner

pure function cross(x, y) result(z)

    real(dp), intent(in) :: x(3), y(3)
    real(dp) :: z(3)

    z = [x(2)*y(3) - x(3)*y(2), x(3)*y(1) - x(1)*y(3), x(1)*y(2) - x(2)*y(1)]

end function cross

end module homma_collision_mod
