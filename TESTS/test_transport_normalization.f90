program test_transport_normalization

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use flux_deviation_mod, only: lorentz_pitch_step
    use gorilla_applets_sub_mod, only: collision_frequency_from_nu_star, &
        & nu_star_from_collision_frequency,radial_metric_ds_dr

    implicit none

    real(dp), parameter :: nu_star = 0.125_dp
    real(dp), parameter :: major_radius = 165.0_dp
    real(dp), parameter :: aiota = -0.4_dp
    real(dp), parameter :: speed = 5.0e7_dp
    real(dp), parameter :: pitch = 0.3_dp
    real(dp), parameter :: eps_collision = 0.02_dp
    real(dp) :: nu_collision,pitch_plus,pitch_minus,expected_second_moment
    real(dp) :: ds_dr_positive,ds_dr_negative

    nu_collision = collision_frequency_from_nu_star(nu_star,major_radius,aiota,speed)
    if(abs(nu_collision-15151.515151515152_dp) > 1.0e-9_dp) then
        error stop 'standard nu-star conversion has the wrong normalization'
    endif
    if(abs(nu_star_from_collision_frequency(nu_collision,major_radius,aiota,speed) &
        & - nu_star) > 1.0e-15_dp) then
        error stop 'nu-star conversion does not round trip'
    endif

    ds_dr_positive = radial_metric_ds_dr(0.49_dp,2.0e4_dp,3.5e7_dp)
    ds_dr_negative = radial_metric_ds_dr(0.49_dp,2.0e4_dp,-3.5e7_dp)
    if(ds_dr_positive <= 0.0_dp .or. ds_dr_positive /= ds_dr_negative) then
        error stop 'radial metric must be positive and invariant to flux handedness'
    endif

    pitch_plus = lorentz_pitch_step(pitch,eps_collision,1.0_dp)
    pitch_minus = lorentz_pitch_step(pitch,eps_collision,-1.0_dp)
    if(abs(0.5_dp*(pitch_plus+pitch_minus)-pitch*(1.0_dp-eps_collision)) &
        & > 1.0e-15_dp) then
        error stop 'Lorentz step does not reproduce the P1 decay oracle'
    endif
    expected_second_moment = pitch**2*(1.0_dp-eps_collision)**2 &
        & +(1.0_dp-pitch**2)*eps_collision
    if(abs(0.5_dp*(pitch_plus**2+pitch_minus**2)-expected_second_moment) &
        & > 1.0e-15_dp) then
        error stop 'Lorentz step has the wrong second moment'
    endif

    print *, 'PASS: standard nu-star and Lorentz moment oracles'

end program test_transport_normalization
