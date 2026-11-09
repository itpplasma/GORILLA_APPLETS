module thermal_force_rng_mod

! Per-particle random numbers for the thermal-force applet: xoshiro256**
! (Blackman & Vigna 2018) seeded by SplitMix64 from a base seed and the
! particle index. Each particle owns its stream, so results do not depend
! on the OpenMP schedule, and a gradient run and its control run can use
! identical numbers. Relies on int64 wrap-around: compile with -fwrapv.

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64

    implicit none

    private

    type, public :: rng_t
        integer(int64) :: s(4) = 0_int64
    end type rng_t

    real(dp), parameter :: two_pi = 6.28318530717958647692528676655900577_dp

    public :: rng_seed, rng_uniform, rng_normal3

contains

subroutine rng_seed(rng, seed, index)

    type(rng_t), intent(out) :: rng
    integer(int64), intent(in) :: seed, index
    integer(int64) :: sm, mixed
    integer :: k

    sm = seed
    mixed = splitmix64(sm)
    sm = ieor(mixed, index)
    do k = 1, 4
        rng%s(k) = splitmix64(sm)
    end do

end subroutine rng_seed

function rng_uniform(rng) result(u)

    ! Uniform on [0, 1) with 53 random bits.

    type(rng_t), intent(inout) :: rng
    real(dp) :: u

    u = real(ishft(next(rng), -11), dp)*2.0_dp**(-53)

end function rng_uniform

subroutine rng_normal3(rng, xi)

    ! Three independent standard normal numbers (Box-Muller).

    type(rng_t), intent(inout) :: rng
    real(dp), intent(out) :: xi(3)
    real(dp) :: r, phi

    r = sqrt(-2.0_dp*log(1.0_dp - rng_uniform(rng)))
    phi = two_pi*rng_uniform(rng)
    xi(1) = r*cos(phi)
    xi(2) = r*sin(phi)
    r = sqrt(-2.0_dp*log(1.0_dp - rng_uniform(rng)))
    phi = two_pi*rng_uniform(rng)
    xi(3) = r*cos(phi)

end subroutine rng_normal3

function next(rng) result(r)

    type(rng_t), intent(inout) :: rng
    integer(int64) :: r, t

    r = ishftc(rng%s(2)*5_int64, 7)*9_int64
    t = ishft(rng%s(2), 17)
    rng%s(3) = ieor(rng%s(3), rng%s(1))
    rng%s(4) = ieor(rng%s(4), rng%s(2))
    rng%s(2) = ieor(rng%s(2), rng%s(3))
    rng%s(1) = ieor(rng%s(1), rng%s(4))
    rng%s(3) = ieor(rng%s(3), t)
    rng%s(4) = ishftc(rng%s(4), 45)

end function next

function splitmix64(x) result(z)

    integer(int64), intent(inout) :: x
    integer(int64) :: z

    x = x - 7046029254386353131_int64            ! + 0x9E3779B97F4A7C15
    z = x
    z = ieor(z, ishft(z, -30))*(-4658895280553007687_int64)   ! 0xBF58476D1CE4E5B9
    z = ieor(z, ishft(z, -27))*(-7723592293110705685_int64)   ! 0x94D049BB133111EB
    z = ieor(z, ishft(z, -31))

end function splitmix64

end module thermal_force_rng_mod
