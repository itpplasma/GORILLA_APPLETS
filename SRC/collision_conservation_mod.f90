module collision_conservation_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: update_background_reservoir

contains

subroutine update_background_reservoir(test_mass, background_mass, &
        marker_to_background_ratio, delta_test_vpar, delta_test_energy, &
        background_vpar, background_temperature)

    real(dp), intent(in) :: test_mass
    real(dp), intent(in) :: background_mass
    real(dp), intent(in) :: marker_to_background_ratio
    real(dp), intent(in) :: delta_test_vpar
    real(dp), intent(in) :: delta_test_energy
    real(dp), intent(inout) :: background_vpar
    real(dp), intent(inout) :: background_temperature

    real(dp) :: background_vpar_old

    if (test_mass <= 0.0_dp) error stop 'test_mass must be positive'
    if (background_mass <= 0.0_dp) error stop 'background_mass must be positive'
    if (marker_to_background_ratio < 0.0_dp) &
        error stop 'marker_to_background_ratio must be non-negative'

    background_vpar_old = background_vpar
    background_vpar = background_vpar_old - marker_to_background_ratio* &
        test_mass*delta_test_vpar/background_mass
    background_temperature = background_temperature - &
        (2.0_dp/3.0_dp)*marker_to_background_ratio*delta_test_energy - &
        (background_mass/3.0_dp)*(background_vpar**2 - &
        background_vpar_old**2)

end subroutine update_background_reservoir

end module collision_conservation_mod
