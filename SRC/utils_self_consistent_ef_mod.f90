module utils_self_consistent_ef_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine read_self_consistent_electric_field_inp_into_type

    use gorilla_applets_types_mod, only: in

    real(dp) :: time_step,energy_eV,n_particles, density
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_boltzmann_energies, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
               boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
               boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, boole_static_ne
    integer :: i_integrator_type, seed_option, n_electric_potential_updates, update_dimension, n_species

    integer :: s_inp_unit

    !Namelist for self consistent electric field input
    NAMELIST /self_consistent_ef_nml/ time_step,energy_eV,n_particles,boole_squared_moments,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, &
    & boole_antithetic_variate,boole_linear_temperature_simulation,i_integrator_type,seed_option, boole_write_vertex_indices, &
    & boole_write_vertex_coordinates, boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
    & boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_electric_potential_updates, update_dimension, &
    & n_species, boole_static_ne

    open(newunit = s_inp_unit, file='self_consistent_ef.inp', status='unknown')
    read(s_inp_unit,nml=self_consistent_ef_nml)
    close(s_inp_unit)

    in%time_step = time_step
    in%energy_eV = energy_eV
    in%n_particles = n_particles
    in%density = density
    in%boole_squared_moments = boole_squared_moments
    in%boole_point_source = boole_point_source
    in%boole_collisions = boole_collisions
    in%boole_precalc_collisions = boole_precalc_collisions
    in%boole_refined_sqrt_g = boole_refined_sqrt_g
    in%boole_boltzmann_energies = boole_boltzmann_energies
    in%boole_linear_density_simulation = boole_linear_density_simulation
    in%boole_antithetic_variate = boole_antithetic_variate
    in%boole_linear_temperature_simulation = boole_linear_temperature_simulation
    in%i_integrator_type = i_integrator_type
    in%seed_option = seed_option
    in%num_particles = int(n_particles)
    in%boole_write_vertex_indices = boole_write_vertex_indices
    in%boole_write_vertex_coordinates = boole_write_vertex_coordinates
    in%boole_write_prism_volumes = boole_write_prism_volumes
    in%boole_write_refined_prism_volumes = boole_write_refined_prism_volumes
    in%boole_write_boltzmann_density = boole_write_boltzmann_density
    in%boole_write_electric_potential = boole_write_electric_potential
    in%boole_write_moments = boole_write_moments
    in%boole_write_fourier_moments = boole_write_fourier_moments
    in%boole_write_exit_data = boole_write_exit_data
    in%boole_write_grid_data = boole_write_grid_data
    in%boole_preserve_energy_and_momentum_during_collisions = boole_preserve_energy_and_momentum_during_collisions
    in%n_electric_potential_updates = n_electric_potential_updates
    in%update_dimension = update_dimension
    in%n_species = n_species
    in%boole_static_ne = boole_static_ne

    print *,'GORILLA_APPLETS: Loaded input data from self_consistent_ef.inp'

end subroutine read_self_consistent_electric_field_inp_into_type

subroutine allocate_electric_potential_type

    use gorilla_applets_types_mod, only: ep, in
    use tetra_grid_mod, only: ntetr, nvert
    use tetra_grid_settings_mod, only: grid_size

    allocate(ep%rho_prism(ntetr/3))
    allocate(ep%rho_flux_layer(grid_size(1)))
    allocate(ep%rho_vert(nvert))
    allocate(ep%phi_elec_from_rho(nvert))
    allocate(ep%average_abs_phi_elec_from_rho(max(1,in%n_electric_potential_updates)))
    allocate(ep%total_tracing_time(max(1,in%n_electric_potential_updates)))

end subroutine allocate_electric_potential_type

subroutine perform_electric_potential_update(i)

    use gorilla_applets_types_mod, only: c, output, grid_t, in, ep, one_d
    use gorilla_settings_mod, only: boole_time_Hamiltonian, coord_system
    use tetra_grid_mod, only: ntetr, nvert
    use tetra_physics_mod, only: make_tetra_physics, phi_elec
    use field_mod, only: ipert

    integer, intent(in) :: i
    real(dp) :: r=0.0_dp !under-relaxation factor
    logical :: boole_print_1d_data = .false.

    if (boole_time_Hamiltonian.eqv..false.) then
        print*, "Error, variable 'boole_time_Hamiltonian' must be set to '.true.' for electric potential update to work"
        stop
    endif

    one_d%boole_print_densities = .true.

    call calc_phi_elec_from_rho(i)

    phi_elec = phi_elec + ep%phi_elec_from_rho

    !giving i_option = 12 as input variable prevents make_tetra_physics from overwriting phi_elec
    call make_tetra_physics(coord_system,ipert,i_option = 12)

    if (i.gt.0) call print_data(i)

end subroutine perform_electric_potential_update

subroutine calc_phi_elec_from_rho(i)

    use gorilla_applets_types_mod, only: in, ep, start, output
    use constants, only: ev2erg, eps, echarge

    integer, intent(in) :: i
    real(dp) :: factor

    ep%rho_prism = 0
    ep%rho_flux_layer = 0
    ep%rho_vert = 0
    
    ep%phi_elec_from_rho = 0
    if (i.gt.0) then
        ep%average_abs_phi_elec_from_rho(i) = 0
        ep%total_tracing_time(i) = 0
    endif
    
    if (in%update_dimension.eq.1) then

        call calc_average_charge_density_per_flux_layer(i)
        call calc_rho_on_vertices

        ! if (i.eq.1) ep%mean_abs_rho_at_first_update = sum(abs(ep%rho_vert))/size(ep%rho_vert)

        ! if (ep%mean_abs_rho_at_first_update.gt.eps**2) then
        !     factor = in%energy_eV*ev2erg/echarge/ep%mean_abs_rho_at_first_update
        ! else
        !     factor = 1.0_dp
        ! endif

        
        factor = in%energy_eV*ev2erg/(echarge**2*in%density)!T/(n_0*e^2)= 4*pi*r_D^2

        ep%phi_elec_from_rho =  ep%rho_vert*factor
        
    endif

end subroutine calc_phi_elec_from_rho

subroutine calc_average_charge_density_per_flux_layer(i)

    use gorilla_applets_types_mod, only:  g, output, ep, in, start, one_d, exit_data
    use tetra_grid_mod, only: nvert, verts_sthetaphi
    use tetra_grid_settings_mod, only: grid_size
    use constants, only: echarge

    integer, intent(in) :: i 
    integer :: ns, j, species, ind_prism, n_species
    real(dp), dimension(:), allocatable :: s_shell_volumes, electron_densities
    real(dp) :: s, electron_density_factor, factor_from_ion_weights

    n_species = in%n_species

    allocate(s_shell_volumes(grid_size(1)),electron_densities(grid_size(1)))
    if (one_d%boole_print_densities) then
        if (.not.allocated(one_d%densities)) allocate(one_d%densities(grid_size(1),in%n_species))
        one_d%densities = 0.0_dp
    endif

    s_shell_volumes = 0.0_dp

    do j = 1,in%num_particles
        do species = 1,2
            ep%total_tracing_time(i) = ep%total_tracing_time(i) + exit_data%t_confined(j,species)
        enddo
    enddo

    !Compute s_shell_volumes and if in%boole_static_ne compute electron_densities
    do ns = 1,grid_size(1)
        s = (verts_sthetaphi(1, grid_size(3)*(ns-1)+1) + verts_sthetaphi(1,grid_size(3)*ns+1))/2.0_dp
        if (in%boole_static_ne) electron_densities(ns) = in%density*(1.0_dp-0.9_dp*s)
        do j = 1, grid_size(2)*grid_size(3)*2
            ind_prism = g%prisms_per_flux_tube(ns,j)
            s_shell_volumes(ns) = s_shell_volumes(ns) + output%prism_volumes(ind_prism)
        enddo
    enddo

    if (in%boole_static_ne) then
        n_species = in%n_species-1
        if (in%boole_linear_density_simulation) then
            electron_density_factor = 1.0d0
        else
            factor_from_ion_weights = sum(start%weight(:,1))/(in%num_particles*in%density*sum(output%prism_volumes(:)))
            electron_density_factor = in%density/(sum(electron_densities*s_shell_volumes)/sum(s_shell_volumes))*&
                                      factor_from_ion_weights
        endif
    endif

    do ns = 1,grid_size(1)
        s = (verts_sthetaphi(1, grid_size(3)*(ns-1)+1) + verts_sthetaphi(1,grid_size(3)*ns+1))/2.0_dp
        do j = 1, grid_size(2)*grid_size(3)*2
            ind_prism = g%prisms_per_flux_tube(ns,j)
            do species = 1,n_species
                ep%rho_flux_layer(ns) = ep%rho_flux_layer(ns) + real(output%prism_moments(1,ind_prism,species))* &
                                        output%prism_volumes(ind_prism)*start%particle_charge(species)
                if (one_d%boole_print_densities) one_d%densities(ns,species) = one_d%densities(ns,species) + &
                                                    real(output%prism_moments(1,ind_prism,species))*output%prism_volumes(ind_prism)
            enddo
            if (in%boole_static_ne) then
                ep%rho_flux_layer(ns) = ep%rho_flux_layer(ns) + electron_density_factor*electron_densities(ns)*&
                                        output%prism_volumes(ind_prism)*(-echarge)
                                        !*ep%total_tracing_time(i)/(in%num_particles*in%time_step)
                if (one_d%boole_print_densities) one_d%densities(ns,in%n_species) = one_d%densities(ns,in%n_species) + &
                                                 electron_density_factor*electron_densities(ns)*output%prism_volumes(ind_prism)
                                                 !*ep%total_tracing_time(i)/(in%num_particles*in%time_step)
            endif
        enddo
        ep%rho_flux_layer(ns) = ep%rho_flux_layer(ns)/s_shell_volumes(ns)
        if (one_d%boole_print_densities) one_d%densities(ns,:) = one_d%densities(ns,:)/s_shell_volumes(ns)
    enddo

end subroutine calc_average_charge_density_per_flux_layer

subroutine calc_rho_on_vertices

    use tetra_grid_mod, only: nvert, verts_rphiz, verts_sthetaphi
    use tetra_grid_settings_mod, only: grid_size
    use gorilla_applets_types_mod, only: g, ep
    use tetra_physics_mod, only: tetra_physics, coord_system, mag_axis_R0, mag_axis_Z0
    use constants, only: pi

    integer :: ns
    real(dp) :: value_to_be_set
    real(dp) :: rho_surface_2, rho_surface_3, rho_surface_end_minus_1, rho_surface_end_minus_2
    real(dp) :: distance_a, distance_b, extrapolation_factor
    real(dp), dimension(grid_size(1)) :: delta_s
    real(dp), dimension(grid_size(1)+1) :: rho_per_flux_surface

    !compute delta s between consecutive flux surfaces
    do ns = 1,grid_size(1)
        if (coord_system.eq.2) then
            delta_s(ns) = (verts_sthetaphi(1,grid_size(3)*ns+1)-verts_sthetaphi(1,grid_size(3)*(ns-1)+1))
        else 
            delta_s(ns) = sqrt((verts_rphiz(1,grid_size(3)*ns+1)-verts_rphiz(1,grid_size(3)*(ns-1)+1))**2 + &
                               (verts_rphiz(3,grid_size(3)*ns+1)-verts_rphiz(3,grid_size(3)*(ns-1)+1))**2)
        endif
    enddo

    !set rho on even flux surfaces
    do ns = 2,grid_size(1),2
        value_to_be_set = 0.5_dp*(ep%rho_flux_layer(ns-1)+ep%rho_flux_layer(ns))
        !remove this line in the future >
        !value_to_be_set = sin(dble(ns-1)/dble(grid_size(1)-1)*pi)
        !value_to_be_set = dble(ns)/dble(grid_size(1))
        call fill_vector_parts_with_value(ep%rho_vert, g%vertices_per_flux_surface(ns,:), value_to_be_set)
        rho_per_flux_surface(ns) = value_to_be_set
    enddo

    !set rho on odd flux surfaces (without borders)
    do ns = 3,grid_size(1)-1,2
        value_to_be_set = (rho_per_flux_surface(ns-1)*delta_s(ns) + rho_per_flux_surface(ns+1)*delta_s(ns-1))/ &
                          (delta_s(ns-1)+delta_s(ns))
        call fill_vector_parts_with_value(ep%rho_vert, g%vertices_per_flux_surface(ns,:), value_to_be_set)
        rho_per_flux_surface(ns) = value_to_be_set
    enddo

    !set rho on first flux surface
    value_to_be_set = rho_per_flux_surface(2) + (rho_per_flux_surface(2)-rho_per_flux_surface(3))*(delta_s(1)/delta_s(2))
    call fill_vector_parts_with_value(ep%rho_vert, g%vertices_per_flux_surface(1,:), value_to_be_set)
    rho_per_flux_surface(1) = value_to_be_set

    if (mod(grid_size(1),2).eq.1) then !there is an odd number of flux surfaces and the last but one flux surface has to be set
        value_to_be_set = rho_per_flux_surface(grid_size(1)-1) + (delta_s(grid_size(1)-1)/delta_s(grid_size(1)-2))* &
                         (rho_per_flux_surface(grid_size(1)-1)-rho_per_flux_surface(grid_size(1)-2))
        call fill_vector_parts_with_value(ep%rho_vert, g%vertices_per_flux_surface(grid_size(1),:), value_to_be_set)
        rho_per_flux_surface(grid_size(1)) = value_to_be_set
    endif

    !set rho on last flux surface
    value_to_be_set = rho_per_flux_surface(grid_size(1)) + (delta_s(grid_size(1))/delta_s(grid_size(1)-1))* &
                     (rho_per_flux_surface(grid_size(1))-rho_per_flux_surface(grid_size(1)-1))
    call fill_vector_parts_with_value(ep%rho_vert, g%vertices_per_flux_surface(grid_size(1)+1,:), value_to_be_set)
    rho_per_flux_surface(grid_size(1)+1) = value_to_be_set

end subroutine calc_rho_on_vertices

subroutine print_data(i)

    use gorilla_applets_types_mod, only: c, output, grid_t, in, ep, exit_data, counter, s, one_d
    use tetra_physics_mod, only: phi_elec
    use tetra_grid_settings_mod, only: grid_size
    use tetra_physics_mod, only: tetra_physics

    integer, intent(in) :: i
    integer :: ep_unit, ed_unit, pe_unit, id_unit, l_unit, ef_unit, s_unit, one_d_unit, j, species, k
    character(len=100) :: filename_ep, filename_ed, i_str, filename_phi_elec, filename_ion_densities, filename_lost, filename_ef, &
                          filename_s, filename_1d

    if (i.eq.1) then
        filename_ef = 'electric_field.dat'
        open(newunit = ef_unit, file = filename_ef)
        do j = 1,grid_size(1)
            write(ef_unit,*) tetra_physics((j-1)*6*grid_size(2)+1)%Er_mod
        enddo
        close(ef_unit)
    endif
    
    filename_lost = 'number_of_lost_particles.dat'
    if (i.eq.1) call unlink(filename_lost)
    open(newunit = l_unit, file = filename_lost,position = 'append')
    write(l_unit,*) i, counter%lost_particles
    close(l_unit)

    write(i_str, '(I0)') i  ! Convert integer to string without leading spaces

    if (one_d%boole_print_densities.eqv..true.) then
        filename_1d = 'one_d_densities' // trim(i_str) // '.dat'
        call unlink(filename_1d)
        open(newunit = one_d_unit, file = filename_1d)
        do j = 1,grid_size(1)
            write(one_d_unit,*) one_d%densities(j,:)
        enddo
        close(one_d_unit)
    endif

    filename_ep = 'rho_per_vertex_during_electric_potential_update_' // trim(i_str) // '.dat'
    call unlink(filename_ep)
    !open(newunit = ep_unit, file = filename_ep)
    !write(ep_unit,'(ES20.10E4)') ep%rho_vert
    !close(ep_unit)

    ! filename_s = 's_statistics_after_electric_potential_update' // trim(i_str) // '.dat'
    ! call unlink(filename_s)
    ! k = size(s%time)
    ! open(newunit = s_unit, file = filename_s)
    ! do j = 1,k
    !     write(s_unit,*) s%time(j), s%delta_s(j), s%delta_s_squared(j), s%check(j)
    ! enddo
    ! close(s_unit)

    filename_phi_elec = 'phi_elec_after_electric_potential_update_' // trim(i_str) // '.dat'
    call unlink(filename_phi_elec)
    open(newunit = pe_unit, file = filename_phi_elec)
    write(pe_unit,'(ES20.10E4)') phi_elec
    close(pe_unit)

    filename_ion_densities = 'ion_densities_after_electric_potential_update_' // trim(i_str) // '.dat'
    call unlink(filename_ion_densities)
    open(newunit = id_unit, file = filename_ion_densities)
    write(id_unit,'(ES20.10E4)') real(output%prism_moments(1,:,1))
    close(id_unit)

    filename_ed = 'exit_data_' // trim(i_str) // '.dat'
    call unlink(filename_ed)
    open(newunit = ed_unit, file = filename_ed)
    do j=1,in%num_particles
        write(ed_unit,*) exit_data%t_confined(j,1), dble(exit_data%integration_step(j,1)), exit_data%x(:,j,1)
    enddo
    ! do j=1,in%num_particles
    !     write(ed_unit,*) exit_data%t_confined(j,2), dble(exit_data%integration_step(j,2)), exit_data%x(:,j,2)
    ! enddo
    close(ed_unit)



    ep%average_abs_phi_elec_from_rho(i) = sum(abs(ep%phi_elec_from_rho))/size(ep%phi_elec_from_rho)

    print*, "electric potential update ", i, " complete."
    print*, "Average abs(Delta Phi) is ", ep%average_abs_phi_elec_from_rho(i)
    print*, "Maximum abs(Delta Phi) is ", maxval(abs(ep%phi_elec_from_rho))
    print*, "Total tracing time is ", ep%total_tracing_time(i)

end subroutine print_data

subroutine associate_flux_labels_with_tetrahedra_and_vertices

    use tetra_grid_settings_mod, only: grid_kind, grid_size
    use gorilla_applets_types_mod, only: g

    integer :: ns, ntheta, nphi, i, ind_prism, j, ind_vert

    ! if ((.not.grid_kind.eq.3).and.(.not.grid_kind.eq.4))  then
    !     print*, 'Error, flux labels can ony be associated with tetrahedra and vertices if grid_kind is either 3 or 4. &
    !              Program is terminated.'
    !     stop
    ! endif

    allocate(g%vertices_per_flux_surface(grid_size(1)+1,grid_size(2)*grid_size(3)))
    allocate(g%prisms_per_flux_tube(grid_size(1),grid_size(2)*grid_size(3)*2))

    !Fill vertices_per_flux_surface
    do ns = 1,grid_size(1)+1
        i = 1
        do nphi = 1,grid_size(2)
            do ntheta = 1,grid_size(3)
                ind_vert = (nphi-1)*(grid_size(1)+1)*grid_size(3) + (ns-1)*grid_size(3) + ntheta
                g%vertices_per_flux_surface(ns,i) = ind_vert
                i = i + 1
            enddo
        enddo
    enddo

    !Fill prisms_per_flux_tube
    do ns = 1,grid_size(1)
        i = 1
        do nphi = 1,grid_size(2)
            do ntheta = 1,grid_size(3)
                do j = 1,2 !count through prisms per hexahedron
                    ind_prism = (nphi-1)*grid_size(1)*grid_size(3)*2 + (ns-1)*grid_size(3)*2 + (ntheta-1)*2 + j
                    g%prisms_per_flux_tube(ns,i) = ind_prism
                    i = i + 1
                enddo
            enddo
        enddo
    enddo

end subroutine associate_flux_labels_with_tetrahedra_and_vertices

subroutine mirror_particles_on_domain_boundaries(x,vpar,n,ind_tetr,iface,z_save,perpinv,ind_tetr_save)

    use supporting_functions_mod, only: vperp_func
    use tetra_grid_settings_mod, only: sfc_s_min, n_field_periods
    use find_tetra_mod, only: find_tetra
    use gorilla_settings_mod, only: poly_order
    use constants, only: pi

    real(dp), dimension(3), intent(inout) :: x, z_save
    real(dp), intent(in) :: perpinv
    real(dp), intent(inout) :: vpar
    integer, intent(in) :: n, ind_tetr_save
    integer, intent(inout) :: ind_tetr, iface
    real(dp) :: vperp
    real(dp), dimension(3) :: x_new
    logical :: boole_diag = .true.

    x_new = (/x(1),-x(2)+2*pi,-x(3)+2*pi/n_field_periods/)
    vpar = -vpar
    vperp = vperp_func(z_save,perpinv,ind_tetr_save)
    call find_tetra(x_new,vpar,vperp,ind_tetr,iface)

    if (.not.(x(1).gt.sfc_s_min)) then
        !if (boole_diag) print*, "particle ", n, " is being pushed across the central annulus at s = ", x(1)
    else
        !if (boole_diag) print*, "particle ", n, " is being mirrored at s = ", x(1)
    endif
    if (ind_tetr.eq.-1) then
        if (boole_diag) print*, "ATTENTION: particle pushing was unsuccessful, vperp "
        if (boole_diag) print*, "x = ", x_new
        if (boole_diag) print*, "vpar, vperp = ", vpar, vperp
    else
        !if (boole_diag) print*, "particle pushing was successful"
        x = x_new
    endif

end subroutine mirror_particles_on_domain_boundaries

subroutine print_errors_for_bad_inputs

    use gorilla_applets_types_mod, only: in
    use tetra_physics_mod, only: coord_system

    if (in%i_integrator_type.eq.2) then
        print*, 'Error: i_integrator_type set to 2, this module only works with i_integrator_type set to 1'
        print*, 'Program terminated'
        stop
    endif

    if (in%boole_refined_sqrt_g.eqv..true.) then
        print*, 'Error: boole_refined_sqrt_g set to .true., but this only works for cylindrical coordinates. This module &
                 works with flux coordinates.'
        print*, 'Program terminated'
        stop
    endif

    if (in%boole_write_refined_prism_volumes.eqv..true.) then
        print*, 'Error: boole_write_refined_prism_volumes set to .true., but this only works for cylindrical coordinates. &
                 This module works with flux coordinates.'
        print*, 'Program terminated'
        stop
    endif

    if (coord_system.eq.1) then
        print*, 'Error: coord_system set to 1, but this module works with flux coordinates.'
        print*, 'Program terminated'
        stop
    endif

end subroutine print_errors_for_bad_inputs

subroutine fill_vector_parts_with_value(vector,indices,set_value)

    real(dp), dimension(:), intent(inout) :: vector
    integer, dimension(:), intent(in) :: indices
    real(dp), intent(in) :: set_value
    integer :: i, numel

    numel = size(indices)

    do i = 1,numel
        vector(indices(i)) = set_value
    enddo        

end subroutine fill_vector_parts_with_value

subroutine treat_particles_that_are_lost_but_should_not_be(z_save_at_x_save,ind_tetr_save,z_save,x_save,x,vpar,vperp, &
                                                                                                        perpinv,ind_tetr,vpar_save)

    use utils_orbit_timestep_mod, only: initialize_constants_of_motion
    use supporting_functions_mod, only: vperp_func

    integer, intent(in)                          :: ind_tetr_save
    real(dp), dimension(3), intent(in)           :: z_save_at_x_save, x_save
    real(dp), dimension(3), intent(inout)        :: x, z_save
    real(dp), intent(inout)                      :: vpar,vperp, perpinv, vpar_save
    integer, intent(inout)                       :: ind_tetr
    real(dp)                                     :: vperp_save
    real(dp)                                     :: v
    integer                                      :: problem_unit

    print*, 'This should not happen.'
    vperp_save = vperp_func(z_save_at_x_save,perpinv,ind_tetr_save)
    vperp      = vperp_func(z_save,          perpinv,ind_tetr_save)

    print*, 'x_save, vpar_save, vperp_save = ', x_save, vpar_save, vperp_save
    print*, 'x, vpar, vperp = ', x, vpar, vperp

    open(newunit = problem_unit, file = 'pushing_problems.dat', position = 'append')
    write(problem_unit,*) 'x_save, vpar_save, vperp_save = ', x_save, vpar_save, vperp_save
    write(problem_unit,*) 'x, vpar, vperp = ', x, vpar, vperp
    close(problem_unit)

    !Continue tracing the particle from the previous tetrahedron crossing:
    !this time half the value of vpar to avoid running into the same problem again
    v = sqrt(vpar_save**2+vperp_save**2)
    vpar = 0.5_dp*vpar_save
    vperp = sqrt(v**2-vpar**2)
    call initialize_constants_of_motion(vperp,z_save_at_x_save,ind_tetr_save,perpinv)
    x = x_save
    ind_tetr = ind_tetr_save
    z_save = z_save_at_x_save

end subroutine treat_particles_that_are_lost_but_should_not_be

end module utils_self_consistent_ef_mod