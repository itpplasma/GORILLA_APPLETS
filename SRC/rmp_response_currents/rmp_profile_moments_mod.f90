module rmp_profile_moments_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: deposit_profile_moments
contains
    pure subroutine deposit_profile_moments(moments, n_bins, bin, phase, current, density)
        complex(dp), intent(inout) :: moments(:)
        integer, intent(in) :: n_bins, bin
        real(dp), intent(in) :: phase
        complex(dp), intent(in) :: current, density
        complex(dp) :: demodulation

        demodulation = exp(cmplx(0.0_dp, -phase, kind=dp))
        moments(bin) = moments(bin) + current * demodulation
        ! Keep the legacy axisymmetric density diagnostic in the second block.
        moments(n_bins + bin) = moments(n_bins + bin) + density
        moments(2 * n_bins + bin) = moments(2 * n_bins + bin) + density * demodulation
    end subroutine deposit_profile_moments
end module rmp_profile_moments_mod
