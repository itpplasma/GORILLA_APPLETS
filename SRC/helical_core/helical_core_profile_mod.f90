module helical_core_profile_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: linear_density_factor, delta_f_density_weight

contains

    pure real(dp) function linear_density_factor(s, s_zero)
        real(dp), intent(in) :: s, s_zero

        linear_density_factor = (s_zero - s)/s_zero
    end function linear_density_factor

    pure real(dp) function delta_f_density_weight(base_weight, velocity_s, &
            s_zero)
        real(dp), intent(in) :: base_weight, velocity_s, s_zero

        delta_f_density_weight = base_weight*velocity_s/s_zero
    end function delta_f_density_weight

end module helical_core_profile_mod
