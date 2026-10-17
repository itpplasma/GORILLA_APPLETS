module utils_orbit_timestep_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine identify_particles_entering_annulus(x,local_counter,boole_lost_inside)

    use tetra_physics_mod, only: tetra_physics
    use gorilla_applets_types_mod, only: counter_t, g

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

subroutine update_local_tetr_moments(local_tetr_moments,ind_tetr,n,optional_quantities,species_in)

    use gorilla_applets_types_mod, only: moment_specs, weights, in
    use gorilla_settings_mod, only: optional_quantities_type

    type(optional_quantities_type), intent(in)   :: optional_quantities
    integer, intent(in)                          :: ind_tetr, n
    integer, intent(in), optional                :: species_in
    complex(dp), dimension(:,:), intent (inout)  :: local_tetr_moments
    integer                                      :: m
    integer                                      :: species = 1

    if (present(species_in)) species = species_in

    do m = 1,moment_specs%n_moments

        select case(moment_specs%moments_selector(m))
            case(1)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                    & weights%w(n,species)*optional_quantities%t_hamiltonian
            case(2)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                    & weights%w(n,species)*optional_quantities%gyrophase
            case(3)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                    & weights%w(n,species)*optional_quantities%vpar_int
            case(4)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                    & weights%w(n,species)*optional_quantities%vpar2_int
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

subroutine calc_particle_weights_and_jperp(n,z_save,vpar,vperp,ind_tetr, species_in)

    !Marker weight and perpendicular invariant, shared by all applets:
    !    w = f * J_x * J_y / (pdf_position * pdf_velocity_space)
    !f is the physical distribution (monoenergetic, or Boltzmann with optional linear density and
    !temperature profiles), J_x and J_y are the Jacobians of position and velocity space, and the pdfs are
    !the sampling densities of the markers, evaluated from start%dist_position, start%dist_lambda and
    !start%dist_energy (these must be initialised, also when the markers are placed deterministically).

    use gorilla_applets_types_mod, only: in, start, weights, g
    use tetra_physics_mod, only: tetra_physics
    use constants, only: ev2erg, pi
    use volume_integrals_and_sqrt_g_mod, only: sqrt_g
    use supporting_functions_mod, only: bmod_func
    use gorilla_settings_mod, only: coord_system
    use marker_distribution_mod, only: evaluate_distribution_3d, evaluate_distribution_1d, pdf_boltzmann, &
                                       linear_profile_factor

    real(dp), intent(in) :: vpar, vperp
    real(dp), dimension(3), intent(in) :: z_save
    integer, intent(in) :: n,ind_tetr
    integer, intent(in), optional :: species_in
    integer :: species
    real(dp) :: r, z, s_value, x_global(3)
    real(dp) :: density, T, m, q, energy_erg, Phi_elec, omega_c, f, J_x, J_y
    real(dp) :: pdf_position, pdf_lambda, pdf_energy_erg, pdf_gyroangle, pdf_velocity_space

    species = 1
    if (present(species_in)) species = species_in

    m = start%particle_mass(species)
    r = z_save(1)
    z = z_save(3)
    x_global = z_save + tetra_physics(ind_tetr)%x1

    !Physical density and temperature, optionally with linear radial profiles
    !(for delta-f, the density gradient enters through adapt_weights_delta_f instead)
    density = in%density
    T = in%energy_eV*ev2erg
    if (in%boole_linear_density_simulation.or.in%boole_linear_temperature_simulation) then
        s_value = calc_radial_coordinate(ind_tetr, z_save)
        if (in%boole_linear_density_simulation.and.(.not.in%boole_delta_f)) density = density*linear_profile_factor(s_value)
        if (in%boole_linear_temperature_simulation) T = T*linear_profile_factor(s_value)
    endif

    !Spatial Jacobian
    if (coord_system.eq.2) then
        J_x = abs(tetra_physics(ind_tetr)%sqg1 + sum(tetra_physics(ind_tetr)%gsqg*z_save))
    elseif (in%boole_refined_sqrt_g) then
        J_x = (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
            & (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))
    else
        J_x = r + tetra_physics(ind_tetr)%x1(1)
    endif

    !Sampling densities of the markers
    pdf_position = evaluate_distribution_3d(start%dist_position, x_global)
    pdf_lambda = evaluate_distribution_1d(start%dist_lambda, start%pitch(n,species))
    pdf_gyroangle = 1.0_dp/(2.0_dp*pi)

    !Physical distribution and velocity-space Jacobian
    if (in%boole_monoenergetic) then
        f = density*pdf_gyroangle/(2.0_dp*start%v0(species))
        J_y = start%v0(species)
        pdf_velocity_space = pdf_lambda*pdf_gyroangle
    else
        q = start%particle_charge(species)
        energy_erg = start%energy(n,species)*ev2erg
        pdf_energy_erg = evaluate_distribution_1d(start%dist_energy, start%energy(n,species))/ev2erg
        Phi_elec = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z_save)
        f = pdf_boltzmann(density,T,m,energy_erg,q,Phi_elec)
        omega_c = bmod_func(z_save,ind_tetr)/start%cm_over_e(species)
        J_y = m**2*omega_c
        pdf_velocity_space = pdf_lambda*pdf_energy_erg*pdf_gyroangle*omega_c*sqrt(m/(2.0_dp*energy_erg))
    endif

    !A point source has a degenerate position distribution (evaluated as 1): the markers represent the
    !particles of the whole sampling box
    if (in%boole_point_source) f = f*((g%amax-g%amin)*(g%cmax-g%cmin)*2*pi)

    weights%w(n,species) = f*J_x*J_y/(pdf_position*pdf_velocity_space)
    start%jperp(n,species) = m*vperp**2*start%cm_over_e(species)/(2*bmod_func(z_save,ind_tetr))*(-1)
    !-1 because of negative gyrophase

end subroutine calc_particle_weights_and_jperp

function calc_radial_coordinate(ind_tetr, z_save) result(s_value)

    !Normalised radial coordinate of a position inside tetrahedron ind_tetr: the flux label s in flux
    !coordinates, the normalised poloidal flux for axisymmetric devices in cylindrical coordinates

    use gorilla_applets_types_mod, only: flux
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind

    integer, intent(in) :: ind_tetr
    real(dp), dimension(3), intent(in) :: z_save
    real(dp) :: s_value

    if (coord_system == 2) then
        s_value = tetra_physics(ind_tetr)%x1(1) + z_save(1)
    else if (grid_kind /= 3) then
        s_value = (tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*z_save)) / flux%poloidal_max
    else
        print*, 'Error in calc_radial_coordinate: Computing radial coordinate from A_phi is only valid for &
                &axisymmetric devices. For stellarators (grid_kind=3), use flux coordinates (coord_system=2).'
        stop
    endif

end function calc_radial_coordinate

subroutine compute_radial_fluxes(ind_tetr_save,ind_tetr,x)

    use gorilla_applets_types_mod, only: output
    use tetra_grid_settings_mod, only: grid_size

    integer, intent(in) :: ind_tetr_save,ind_tetr
    integer, dimension(2) :: ind_r, ind_in_first_phi_slice
    real(dp), dimension(3)        :: x
    integer :: change, boundary_num

    ind_in_first_phi_slice = mod((/ind_tetr_save,ind_tetr/) - 1,grid_size(1)*grid_size(3)*6) + 1
    ind_r = int((ind_in_first_phi_slice - 1)/(grid_size(3)*6)) + 1
    
    if (.not.(ind_r(1).eq.ind_r(2))) then
        if (.not.(ind_tetr.eq.-1)) then
            change = ind_r(2) - ind_r(1)
            boundary_num = max(ind_r(2), ind_r(1))
            !$omp critical
            output%radial_flux(boundary_num) = output%radial_flux(boundary_num) + change
            !$omp end critical
        else
            !$omp critical
            if (ind_r(1).eq.1) then
                output%radial_flux(ind_r(1)) = output%radial_flux(ind_r(1)) - 1
            else
                output%radial_flux(ind_r(1)+1) = output%radial_flux(ind_r(1)+1) + 1
            endif
            !$omp end critical
        endif
    endif

end subroutine compute_radial_fluxes

end module utils_orbit_timestep_mod