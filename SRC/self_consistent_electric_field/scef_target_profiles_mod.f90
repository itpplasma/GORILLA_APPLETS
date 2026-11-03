module scef_target_profiles_mod

    !Target plasma profiles of the self-consistent electric field (SCEF) benchmark, as functions of
    !rad = r/a = sqrt(s). They define the collision background and the targets the source tuning
    !aims for. Leaf module (no dependencies), so general modules can use it without depending on SCEF.

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

pure function target_density(rad) result(n)

    !n(rad) = 0.37e14 * (0.74 + 0.26*(1 - rad^2.5)^1.5 - 0.06*(1 - exp(-4*rad^2)))  [cm^-3]

    real(dp), intent(in) :: rad
    real(dp) :: n

    n = 0.37d14 * ( 0.74_dp                                     &
                  + 0.26_dp*(1.0_dp - rad**2.5_dp)**1.5_dp      &
                  - 0.06_dp*(1.0_dp - exp(-4.0_dp*rad**2)) )

end function target_density

pure function target_ion_temperature(rad) result(T_eV)

    !T_i(rad) = 1770 * (0.06 + 0.94*(1 - rad^2.9)^2)  [eV]

    real(dp), intent(in) :: rad
    real(dp) :: T_eV

    T_eV = 1770.0_dp*(0.06_dp + 0.94_dp*(1.0_dp - rad**2.9_dp)**2)

end function target_ion_temperature

pure function target_electron_temperature(rad) result(T_eV)

    !T_e(rad) = T_i(rad) + 3190 * (1 - rad^2.1)^5.4  [eV]

    real(dp), intent(in) :: rad
    real(dp) :: T_eV

    T_eV = target_ion_temperature(rad) + 3190.0_dp*(1.0_dp - rad**2.1_dp)**5.4_dp

end function target_electron_temperature

pure function target_temperature(rad, species) result(T_eV)

    !Target temperature of SCEF species 1 (ions) or 2 (electrons)  [eV]

    real(dp), intent(in) :: rad
    integer, intent(in) :: species
    real(dp) :: T_eV

    if (species.eq.1) then
        T_eV = target_ion_temperature(rad)
    else
        T_eV = target_electron_temperature(rad)
    endif

end function target_temperature

pure function rad_from_s(s) result(rad)

    !r/a from the toroidal flux label s, clipped to [0,1] so the profiles stay defined
    !for points slightly outside the grid

    real(dp), intent(in) :: s
    real(dp) :: rad

    rad = sqrt(min(max(s, 0.0_dp), 1.0_dp))

end function rad_from_s

end module scef_target_profiles_mod
