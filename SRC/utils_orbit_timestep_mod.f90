module utils_orbit_timestep_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    !Tetrahedra with moment contributions from the current marker, per thread. Only active between
    !start_ and stop_touched_tetr_tracking; then the moments of a marker are transferred and reset
    !for these tetrahedra only instead of the whole grid.
    integer, dimension(:), allocatable :: touched_tetr
    logical, dimension(:), allocatable :: is_touched
    integer :: n_touched = 0
    !$omp threadprivate(touched_tetr, is_touched, n_touched)

contains

subroutine start_touched_tetr_tracking(ntetr)

    integer, intent(in) :: ntetr

    if (allocated(touched_tetr)) deallocate(touched_tetr)
    if (allocated(is_touched)) deallocate(is_touched)
    allocate(touched_tetr(ntetr), is_touched(ntetr))
    is_touched = .false.
    n_touched = 0

end subroutine start_touched_tetr_tracking

subroutine stop_touched_tetr_tracking

    if (allocated(touched_tetr)) deallocate(touched_tetr)
    if (allocated(is_touched)) deallocate(is_touched)
    n_touched = 0

end subroutine stop_touched_tetr_tracking

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

subroutine place_particle_on_grid(x,vpar,vperp,ind_tetr,iface,z_save)

    !Find the tetrahedron of a particle that has no grid position yet
    !ind_tetr = -1 ... particle lies outside the grid, z_save is not set

    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use find_tetra_mod, only: find_tetra
    use tetra_physics_mod, only: tetra_physics

    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(in)                  :: vpar, vperp
    integer, intent(out)                  :: ind_tetr, iface
    real(dp), dimension(3), intent(inout) :: z_save

    call check_coordinate_domain(x) !Check coordinate domain (optionally perform modulo operation)
    call find_tetra(x,vpar,vperp,ind_tetr,iface) !Find tetrahedron index and face index for position x
    if (ind_tetr.ne.-1) z_save = x-tetra_physics(ind_tetr)%x1

end subroutine place_particle_on_grid

subroutine push_within_tetrahedron(ind_tetr,iface,x,vpar,z_save,t_push,t_pass,boole_t_finished,iper_phi,optional_quantities, &
                                 & boole_pusher_failed)

    !Push the particle with the pusher selected by ipusher for at most t_push, until it leaves the tetrahedron
    !boole_pusher_failed = .true. ... the pusher could not push the particle and removed it (ind_tetr = -1)

    use pusher_tetra_rk_mod, only: pusher_tetra_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type

    integer, intent(inout)                        :: ind_tetr, iface
    real(dp), dimension(3), intent(inout)         :: x
    real(dp), intent(inout)                       :: vpar
    real(dp), dimension(3), intent(inout)         :: z_save
    real(dp), intent(in)                          :: t_push
    real(dp), intent(out)                         :: t_pass
    logical, intent(out)                          :: boole_t_finished
    integer, intent(out)                          :: iper_phi
    type(optional_quantities_type), intent(inout) :: optional_quantities
    logical, intent(out)                          :: boole_pusher_failed

    select case(ipusher)
        case(1)
            call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t_push,t_pass,boole_t_finished,iper_phi,boole_pusher_failed)
        case(2)
            call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_push,&
                                    & t_pass,boole_t_finished,iper_phi,optional_quantities,boole_pusher_failed)
    end select

end subroutine push_within_tetrahedron

subroutine push_across_axis_hole(x,vpar,vperp,z_save,perpinv,ind_tetr_save,ind_tetr,iface,local_counter,boole_lost_inside)

    !For a particle that has left the grid: if it entered the hole around the magnetic axis
    !(boole_lost_inside = .true.), reflect it through the axis
    !ind_tetr = -1 afterwards ... the particle is outside the grid (reflection failed or not in the hole),
    !                            x is unchanged

    use find_tetra_mod, only: find_tetra
    use supporting_functions_mod, only: vperp_func
    use gorilla_applets_types_mod, only: counter_t, g

    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(in)                  :: vpar, perpinv
    real(dp), intent(inout)               :: vperp
    real(dp), dimension(3), intent(in)    :: z_save
    integer, intent(in)                   :: ind_tetr_save
    integer, intent(inout)                :: ind_tetr, iface
    type(counter_t), intent(inout)        :: local_counter
    logical, intent(out)                  :: boole_lost_inside
    real(dp), dimension(3)                :: x_new

    call identify_particles_entering_annulus(x,local_counter,boole_lost_inside)
    if (.not.boole_lost_inside) return

    x_new = 3*(/g%raxis,x(2),g%zaxis/) - 2*x
    vperp = vperp_func(z_save,perpinv,ind_tetr_save)
    call find_tetra(x_new,vpar,vperp,ind_tetr,iface)
    if (ind_tetr.eq.-1) then
        print*, "ATTENTION: particle pushing across the hole surrounding the magnetic axis was unsuccessful"
    else
        x = x_new
    endif

end subroutine push_across_axis_hole

subroutine reinject_towards_axis(x,vpar,vperp,z_save,perpinv,ind_tetr_save,ind_tetr,iface)

    !For a particle that has left the grid at the outer boundary: move it (R and Z, phi unchanged)
    !by a random fraction of 0 to 1% of its distance towards the magnetic axis
    !ind_tetr = -1 afterwards ... the new position is outside the grid as well, x is unchanged

    use find_tetra_mod, only: find_tetra
    use supporting_functions_mod, only: vperp_func
    use gorilla_applets_types_mod, only: g

    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(in)                  :: vpar, perpinv
    real(dp), intent(inout)               :: vperp
    real(dp), dimension(3), intent(in)    :: z_save
    integer, intent(in)                   :: ind_tetr_save
    integer, intent(inout)                :: ind_tetr, iface
    real(dp), dimension(3)                :: x_new
    real(dp)                              :: rand_frac

    call random_number(rand_frac)
    rand_frac = 0.01_dp * rand_frac
    x_new(1) = x(1) + rand_frac * (g%raxis - x(1))
    x_new(2) = x(2)
    x_new(3) = x(3) + rand_frac * (g%zaxis - x(3))
    vperp = vperp_func(z_save,perpinv,ind_tetr_save)
    call find_tetra(x_new,vpar,vperp,ind_tetr,iface)
    if (ind_tetr.ne.-1) x = x_new

end subroutine reinject_towards_axis

subroutine update_local_tetr_moments(local_tetr_moments,ind_tetr,n,optional_quantities,species_in)

    use gorilla_applets_types_mod, only: moment_specs, weights, in
    use gorilla_settings_mod, only: optional_quantities_type

    type(optional_quantities_type), intent(in)   :: optional_quantities
    integer, intent(in)                          :: ind_tetr, n
    integer, intent(in), optional                :: species_in
    complex(dp), dimension(:,:), intent (inout)  :: local_tetr_moments
    integer                                      :: m
    integer                                      :: species

    species = 1
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

    if (allocated(is_touched)) then
        if (.not.is_touched(ind_tetr)) then
            is_touched(ind_tetr) = .true.
            n_touched = n_touched + 1
            touched_tetr(n_touched) = ind_tetr
        endif
    endif

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

subroutine calc_particle_weights_and_jperp(n,z_save,vpar,vperp,ind_tetr, species_in, boole_include_potential_in)

    !Marker weight and perpendicular invariant, shared by all applets:
    !    w = f * J_x * J_y / (pdf_position * pdf_velocity_space)
    !f is the physical distribution (monoenergetic, or Boltzmann with optional linear density and
    !temperature profiles), J_x and J_y are the Jacobians of position and velocity space, and the pdfs are
    !the sampling densities of the markers, evaluated from start%dist_position, start%dist_lambda and
    !start%dist_energy (these must be initialised, also when the markers are placed deterministically).
    !With boole_include_potential_in = .false. the electric potential is dropped from the Boltzmann f
    !(used e.g. by SCEF, where the source normalisation would be shifted by Phi).

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
    logical, intent(in), optional :: boole_include_potential_in
    integer :: species
    logical :: boole_include_potential
    real(dp) :: r, z, s_value, x_global(3)
    real(dp) :: density, T, m, q, energy_erg, Phi_elec, omega_c, f, J_x, J_y
    real(dp) :: pdf_position, pdf_lambda, pdf_energy_erg, pdf_gyroangle, pdf_velocity_space

    species = 1
    if (present(species_in)) species = species_in
    boole_include_potential = .true.
    if (present(boole_include_potential_in)) boole_include_potential = boole_include_potential_in

    m = start%particle_mass(species)
    r = z_save(1)
    z = z_save(3)
    x_global = z_save + tetra_physics(ind_tetr)%x1

    !Physical density and temperature, optionally with linear radial profiles
    !(for delta-f, the density gradient enters through adapt_weights_delta_f instead).
    !start%density is per marker: in%density in the base applets, the particle source value in SCEF.
    density = start%density(n,species)
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
        if (boole_include_potential) then
            Phi_elec = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z_save)
        else
            Phi_elec = 0.0_dp
        endif
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

subroutine calc_particle_jperp(n, z_save, vperp, ind_tetr, species_in)

    !Perpendicular invariant J_perp only, for callers that assign the marker weight themselves (e.g.
    !the diffusion calibration, which uses the uniform weight 1/n_particles).

    use gorilla_applets_types_mod, only: start
    use supporting_functions_mod, only: bmod_func

    real(dp), intent(in) :: vperp
    real(dp), dimension(3), intent(in) :: z_save
    integer, intent(in) :: n, ind_tetr
    integer, intent(in), optional :: species_in
    integer :: species

    species = 1
    if (present(species_in)) species = species_in

    start%jperp(n,species) = -start%particle_mass(species)*vperp**2*start%cm_over_e(species)/(2.0_dp*bmod_func(z_save,ind_tetr))
    !sign from the negative gyrophase convention

end subroutine calc_particle_jperp

function calc_radial_coordinate(ind_tetr, z_save) result(s_value)

    !Normalised radial coordinate of a position inside tetrahedron ind_tetr: the flux label s in flux
    !coordinates, the normalised poloidal flux for axisymmetric devices in cylindrical coordinates

    use utils_data_pre_and_post_processing_mod, only: normalised_poloidal_flux
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind

    integer, intent(in) :: ind_tetr
    real(dp), dimension(3), intent(in) :: z_save
    real(dp) :: s_value

    if (coord_system == 2) then
        s_value = tetra_physics(ind_tetr)%x1(1) + z_save(1)
    else if (grid_kind /= 3) then
        s_value = normalised_poloidal_flux(tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*z_save))
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
            !$omp atomic update
            output%radial_flux(boundary_num) = output%radial_flux(boundary_num) + change
        else
            if (ind_r(1).eq.1) then
                !$omp atomic update
                output%radial_flux(ind_r(1)) = output%radial_flux(ind_r(1)) - 1
            else
                !$omp atomic update
                output%radial_flux(ind_r(1)+1) = output%radial_flux(ind_r(1)+1) + 1
            endif
        endif
    endif

end subroutine compute_radial_fluxes

end module utils_orbit_timestep_mod