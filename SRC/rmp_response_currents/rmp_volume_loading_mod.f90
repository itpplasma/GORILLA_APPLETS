module rmp_volume_loading_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: volume_candidate, analytic_window_volume
contains
    pure subroutine volume_candidate(r0, a, s_lo, s_hi, u, x, acceptance)
        real(dp), intent(in) :: r0, a, s_lo, s_hi, u(3)
        real(dp), intent(out) :: x(3), acceptance
        real(dp) :: edge, rho2_lo, rho2_hi, rho, theta, pi

        pi = acos(-1.0_dp)
        edge = a*a/(r0 + sqrt(r0*r0 - a*a))
        rho2_lo = s_lo*edge*(2.0_dp*r0 - s_lo*edge)
        rho2_hi = s_hi*edge*(2.0_dp*r0 - s_hi*edge)
        rho = sqrt(rho2_lo + u(1)*(rho2_hi - rho2_lo))
        theta = 2.0_dp*pi*u(2)
        x(1) = r0 + rho*cos(theta)
        x(2) = 2.0_dp*pi*u(3)
        x(3) = rho*sin(theta)
        ! Uniform cross-section proposals become uniform dV = R dR dZ dphi.
        acceptance = x(1)/(r0 + sqrt(rho2_hi))
    end subroutine volume_candidate

    pure real(dp) function analytic_window_volume(r0, a, s_lo, s_hi) result(v)
        real(dp), intent(in) :: r0, a, s_lo, s_hi
        real(dp) :: edge, rho2_lo, rho2_hi, pi

        pi = acos(-1.0_dp)
        edge = a*a/(r0 + sqrt(r0*r0 - a*a))
        rho2_lo = s_lo*edge*(2.0_dp*r0 - s_lo*edge)
        rho2_hi = s_hi*edge*(2.0_dp*r0 - s_hi*edge)
        v = 2.0_dp*pi*pi*r0*(rho2_hi - rho2_lo)
    end function analytic_window_volume
end module rmp_volume_loading_mod
