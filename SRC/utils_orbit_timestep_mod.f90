module utils_orbit_timestep_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine identify_particles_entering_annulus(x,local_counter,boole_lost_inside)

    use tetra_physics_mod, only: tetra_physics
    use boltzmann_types_mod, only: counter_t, pflux, g

    real(dp), dimension(3), intent(in) :: x
    type(counter_t), intent(inout) :: local_counter
    logical, intent(out)   :: boole_lost_inside
    real(dp) :: distance_from_magnetic_axis

    boole_lost_inside = .false.
    distance_from_magnetic_axis = sqrt((x(1)-g%raxis)**2 + (x(3)-g%zaxis)**2)
    if (distance_from_magnetic_axis.lt.g%dist_from_o_point_within_grid) then
        local_counter%lost_inside = local_counter%lost_inside + 1
        boole_lost_inside = .true.
    endif

end subroutine identify_particles_entering_annulus

subroutine update_local_tetr_moments(local_tetr_moments,ind_tetr,n,optional_quantities)

    use boltzmann_types_mod, only: moment_specs, start
    use gorilla_settings_mod, only: optional_quantities_type

    type(optional_quantities_type), intent(in)   :: optional_quantities
    integer, intent(in)                          :: ind_tetr, n
    complex(dp), dimension(:,:), intent (inout)  :: local_tetr_moments
    integer                                      :: m

    do m = 1,moment_specs%n_moments
        select case(moment_specs%moments_selector(m))
            case(1)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                                & start%weight(n)*optional_quantities%t_hamiltonian
            case(2)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                                & start%weight(n)*optional_quantities%gyrophase*start%jperp(n)
            case(3)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                                & start%weight(n)*optional_quantities%vpar_int
            case(4)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                                & start%weight(n)*optional_quantities%vpar2_int
        end select
    enddo

end subroutine update_local_tetr_moments

subroutine initialize_constants_of_motion(vperp,z_save,ind_tetr,perpinv)

    use pusher_tetra_rk_mod, only: initialize_const_motion_rk
    use pusher_tetra_poly_mod, only: initialize_const_motion_poly
    use gorilla_settings_mod, only: ipusher
    use supporting_functions_mod, only: bmod_func

    real(dp), intent(in)               :: vperp
    real(dp), dimension(3), intent(in) :: z_save
    integer, intent(in)                :: ind_tetr
    real(dp), intent(out)              :: perpinv

    perpinv=-0.5_dp*vperp**2/bmod_func(z_save,ind_tetr)
             
    !Initialize constants of motion in particle-private module
    select case(ipusher)
        case(1)
            call initialize_const_motion_rk(perpinv,perpinv**2)
        case(2)
            call initialize_const_motion_poly(perpinv,perpinv**2)
    end select

end subroutine initialize_constants_of_motion

subroutine calc_particle_weights_and_jperp(n,z_save,vpar,vperp,ind_tetr)

    use boltzmann_types_mod, only: in, pflux, start
    use tetra_physics_mod, only: tetra_physics,particle_mass,particle_charge,cm_over_e
    use constants, only: ev2erg
    use volume_integrals_and_sqrt_g_mod, only: sqrt_g
    use supporting_functions_mod, only: bmod_func

    real(dp), intent(in) :: vpar, vperp
    real(dp), dimension(3), intent(in) :: z_save
    integer, intent(in) :: n,ind_tetr
    real(dp) :: local_poloidal_flux, phi_elec_func, temperature
    real(dp) :: r, phi, z

    r = z_save(1)
    phi = z_save(2)
    z = z_save(3)

    if (.not.in%boole_refined_sqrt_g) start%weight = start%weight*(r + tetra_physics(ind_tetr)%x1(1))
    if (in%boole_refined_sqrt_g) then
        start%weight(n) = start%weight(n)* (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
                                        &  (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))               
    endif

    if (in%boole_linear_density_simulation.or.in%boole_linear_temperature_simulation) then
        local_poloidal_flux = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*z_save)
    endif
    if (in%boole_linear_density_simulation) then
        start%weight(n) = start%weight(n)*(pflux%max*1.1_dp-local_poloidal_flux)/(pflux%max*1.1_dp)
    endif

    if (in%boole_boltzmann_energies) then
        !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts have been added before)
        phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z_save)
        if (.not. in%boole_linear_temperature_simulation) then
            start%weight(n) =start%weight(n)*sqrt(start%energy(n)*ev2erg)/(in%energy_eV*ev2erg)**1.5_dp* &
                        & exp(-(start%energy(n)*ev2erg+particle_charge*phi_elec_func)/(in%energy_eV*ev2erg))
        else
            temperature = in%energy_eV*ev2erg*(pflux%max*1.1_dp-local_poloidal_flux)/(pflux%max*1.1_dp)
            start%weight(n) = start%weight(n)*sqrt(start%energy(n)*ev2erg)/temperature**1.5_dp* &
            & exp(-(start%energy(n)*ev2erg+particle_charge*phi_elec_func)/temperature)
        endif
    endif

    start%jperp(n) = particle_mass*vperp**2*cm_over_e/(2*bmod_func(z_save,ind_tetr))*(-1) !-1 because of negative gyrophase

    ! print*, "r sqrt(g) = ", r + tetra_physics(ind_tetr)%x1(1), (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
    ! &  (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))  
    ! if (in%boole_linear_density_simulation) then
    !     print*, "linear density factor = ", (pflux%max*1.1_dp-local_poloidal_flux)/(pflux%max*1.1_dp)
    ! endif
    ! if (in%boole_boltzmann_energies) then
    !     if (.not. in%boole_linear_temperature_simulation) then
    !         print*, "boltzmann energy factor without temperature gradient = ", &
    !         & sqrt(start%energy(n)*ev2erg)/(in%energy_eV*ev2erg)**1.5_dp* &
    !         & exp(-(start%energy(n)*ev2erg+particle_charge*phi_elec_func)/(in%energy_eV*ev2erg))
    !     else
    !         print*, "boltzmann energy factor with temperature gradient = ",sqrt(start%energy(n)*ev2erg)/temperature**1.5_dp* &
    !         & exp(-(start%energy(n)*ev2erg+particle_charge*phi_elec_func)/temperature)
    !     endif
    ! endif

end subroutine calc_particle_weights_and_jperp

end module utils_orbit_timestep_mod