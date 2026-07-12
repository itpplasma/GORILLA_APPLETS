module helical_core_fading_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use constants, only: pi

    implicit none
    private

    public :: delta_f_fade_factor

contains

    pure real(dp) function delta_f_fade_factor(time, total_time, fraction)
        real(dp), intent(in) :: time, total_time, fraction

        real(dp) :: fade_start

        fade_start = total_time*(1.0_dp - fraction)
        if (time <= fade_start) then
            delta_f_fade_factor = 1.0_dp
        else if (time >= total_time) then
            delta_f_fade_factor = 0.0_dp
        else
            delta_f_fade_factor = 0.5_dp*(1.0_dp &
                + cos(pi*(time - fade_start)/(total_time*fraction)))
        end if
    end function delta_f_fade_factor

end module helical_core_fading_mod
