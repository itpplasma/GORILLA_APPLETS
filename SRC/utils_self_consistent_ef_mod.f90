module utils_self_consistent_ef_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

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

    use gorilla_applets_types_mod, only: c, output, grid_t, in, ep
    use gorilla_settings_mod, only: boole_time_Hamiltonian, coord_system
    use tetra_grid_mod, only: ntetr, nvert
    use tetra_physics_mod, only: make_tetra_physics, phi_elec
    use field_mod, only: ipert

    integer, intent(in) :: i
    real(dp) :: r=0.0_dp !under-relaxation factor

    if (boole_time_Hamiltonian.eqv..false.) then
        print*, "Error, variable 'boole_time_Hamiltonian' must be set to '.true.' for electric potential update to work"
        stop
    endif

    call calc_phi_elec_from_rho(i)

    phi_elec = phi_elec + ep%phi_elec_from_rho

    !giving i_option = 12 as input variable prevents make_tetra_physics from overwriting phi_elec
    call make_tetra_physics(coord_system,ipert,i_option = 12)

    if (i.gt.0) call print_data(i)

end subroutine perform_electric_potential_update

subroutine calc_phi_elec_from_rho(i)

    use gorilla_applets_types_mod, only: in, ep
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

        call calc_average_charge_density_per_flux_layer
        call calc_rho_on_vertices

        ! if (i.eq.1) ep%mean_abs_rho_at_first_update = sum(abs(ep%rho_vert))/size(ep%rho_vert)

        ! if (ep%mean_abs_rho_at_first_update.gt.eps**2) then
        !     factor = in%energy_eV*ev2erg/echarge/ep%mean_abs_rho_at_first_update
        ! else
        !     factor = 1.0_dp
        ! endif

        factor = 21.8905_dp
        factor = in%energy_eV*ev2erg/(echarge**2*in%density)!T/(n_0*e^2)= 4*pi*r_D^2

        ep%phi_elec_from_rho =  ep%rho_vert*factor
        
    endif

end subroutine calc_phi_elec_from_rho

subroutine calc_average_charge_density_per_flux_layer

    use gorilla_applets_types_mod, only:  g, output, ep, in, start
    use tetra_grid_mod, only: nvert, verts_sthetaphi
    use tetra_grid_settings_mod, only: grid_size
    use constants, only: echarge

    integer :: ns, j, species, ind_prism, n_species
    real(dp), dimension(:), allocatable :: s_shell_volumes, electron_densities
    real(dp) :: s, electron_density_factor

    n_species = in%n_species

    allocate(s_shell_volumes(grid_size(1)),electron_densities(grid_size(1)))
    s_shell_volumes = 0.0_dp

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
            electron_density_factor = in%density/(sum(electron_densities*s_shell_volumes)/sum(s_shell_volumes))
        endif
    endif

    do ns = 1,grid_size(1)
        s = (verts_sthetaphi(1, grid_size(3)*(ns-1)+1) + verts_sthetaphi(1,grid_size(3)*ns+1))/2.0_dp
        do j = 1, grid_size(2)*grid_size(3)*2
            ind_prism = g%prisms_per_flux_tube(ns,j)
            do species = 1,n_species
                ep%rho_flux_layer(ns) = ep%rho_flux_layer(ns) + real(output%prism_moments(1,ind_prism,species))* &
                                        output%prism_volumes(ind_prism)*start%particle_charge(species)
            enddo
            if (in%boole_static_ne) ep%rho_flux_layer(ns) = ep%rho_flux_layer(ns) + electron_density_factor*electron_densities(ns)*&
                                                            output%prism_volumes(ind_prism)*(-echarge)
        enddo
        ep%rho_flux_layer(ns) = ep%rho_flux_layer(ns)/s_shell_volumes(ns)
    enddo

    print*, 'electron_density_factor = ', electron_density_factor

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

    contains

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

end subroutine calc_rho_on_vertices

subroutine print_data(i)

    use gorilla_applets_types_mod, only: c, output, grid_t, in, ep, exit_data, counter, s
    use tetra_physics_mod, only: phi_elec
    use tetra_grid_settings_mod, only: grid_size
    use tetra_physics_mod, only: tetra_physics

    integer, intent(in) :: i
    integer :: ep_unit, ed_unit, pe_unit, id_unit, l_unit, ef_unit, s_unit, j, species, k
    character(len=100) :: filename_ep, filename_ed, i_str, filename_phi_elec, filename_ion_densities, filename_lost, filename_ef, &
                          filename_s

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

    ! filename_ed = 'exit_data_' // trim(i_str) // '.dat'
    ! call unlink(filename_ed)
    ! open(newunit = ed_unit, file = filename_ed)
    ! do j=1,in%num_particles
    !     write(ed_unit,*) exit_data%t_confined(j,1), dble(exit_data%integration_step(j,1)), exit_data%x
    ! enddo
    ! do j=1,in%num_particles
    !     write(ed_unit,*) exit_data%t_confined(j,2), dble(exit_data%integration_step(j,2)), exit_data%x
    ! enddo
    ! close(ed_unit)

    do j = 1,in%num_particles
        do species = 1,2
            ep%total_tracing_time(i) = ep%total_tracing_time(i) + exit_data%t_confined(j,species)
        enddo
    enddo

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

end module utils_self_consistent_ef_mod