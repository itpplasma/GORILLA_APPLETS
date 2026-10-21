module rmp_maxwellian_loading_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: maxwellian_from_uniforms
contains
    pure subroutine maxwellian_from_uniforms(u, energy_over_t, pitch)
        real(dp), intent(in) :: u(4)
        real(dp), intent(out) :: energy_over_t, pitch
        real(dp) :: z(3), radius, phase, norm2, pi

        pi = acos(-1.0_dp)
        ! For random_number, 1-u is strictly positive, including u=0.
        radius = sqrt(-2.0_dp*log(1.0_dp-u(1)))
        phase = 2.0_dp*pi*u(2)
        z(1) = radius*cos(phase)
        z(2) = radius*sin(phase)
        z(3) = sqrt(-2.0_dp*log(1.0_dp-u(3)))*cos(2.0_dp*pi*u(4))
        norm2 = sum(z*z)
        energy_over_t = 0.5_dp*norm2
        pitch = 0.0_dp
        if (norm2 > 0.0_dp) pitch = z(3)/sqrt(norm2)
    end subroutine maxwellian_from_uniforms
end module rmp_maxwellian_loading_mod
