module utils_data_pre_and_post_processing_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains 

subroutine set_seed_for_random_numbers

    integer,dimension(:), allocatable   :: seed
    integer                             :: seed_inp_unit
    integer                             :: n
    
    open(newunit = seed_inp_unit, file='seed.inp', status='old',action = 'read')
    read(seed_inp_unit,*) n
    allocate(seed(n))
    read(seed_inp_unit,*) seed
    close(seed_inp_unit)
    CALL RANDOM_SEED (PUT=seed)
    deallocate(seed)

end subroutine set_seed_for_random_numbers

subroutine get_ipert()

    use field_mod, only: ipert
    
    integer :: ipert_unit
    
    open(newunit = ipert_unit, file='field_divB0.inp')
    read(ipert_unit,*) ipert        ! 0=eq only, 1=vac, 2=vac+plas no derivatives,
    close(ipert_unit)
    
end subroutine get_ipert

subroutine set_moment_specifications

    use gorilla_settings_mod, only: boole_array_optional_quantities
    use tetra_grid_settings_mod, only: grid_size
    use tetra_grid_mod, only: ntetr
    use gorilla_applets_types_mod, only: moment_specs, in
    
    integer :: i, n_prisms
    
    n_prisms = ntetr/3
    
    moment_specs%boole_squared_moments = in%boole_squared_moments
    moment_specs%n_triangles = n_prisms/grid_size(2)
    moment_specs%n_fourier_modes = 5
    moment_specs%n_moments = 0
    moment_specs%moments_selector = 0
    do i = 1,size(boole_array_optional_quantities)
    if (boole_array_optional_quantities(i).eqv..true.) then
    moment_specs%n_moments = moment_specs%n_moments + 1
    moment_specs%moments_selector(moment_specs%n_moments) = i
    endif
    enddo
    
end subroutine set_moment_specifications

subroutine initialise_output

    use tetra_grid_mod, only: ntetr
    use gorilla_applets_types_mod, only: moment_specs, output
    use tetra_grid_settings_mod, only: grid_size
    
    integer :: n_prisms
    
    n_prisms = ntetr/3
    
    allocate(output%prism_volumes(n_prisms))
    allocate(output%refined_prism_volumes(n_prisms))
    allocate(output%electric_potential(n_prisms))
    allocate(output%boltzmann_density(n_prisms))
    allocate(output%radial_flux(grid_size(1)+1))
    allocate(output%tetr_moments(moment_specs%n_moments,ntetr))
    allocate(output%prism_moments(moment_specs%n_moments,n_prisms))
    if (moment_specs%boole_squared_moments) allocate(output%prism_moments_squared(moment_specs%n_moments,n_prisms))
    allocate(output%moments_in_frequency_space(moment_specs%n_moments,moment_specs%n_triangles,moment_specs%n_fourier_modes))
    
    output%prism_volumes = 0
    output%refined_prism_volumes = 0
    output%electric_potential = 0
    output%boltzmann_density = 0
    output%radial_flux = 0
    output%tetr_moments = 0
    output%prism_moments = 0
    if (moment_specs%boole_squared_moments) output%prism_moments_squared = 0
    output%moments_in_frequency_space = 0
    
end subroutine initialise_output

subroutine calc_starting_conditions(verts)

    use gorilla_applets_types_mod, only: in
    
    real(dp), dimension(:,:), allocatable, intent(out)     :: verts
    real(dp), dimension(:,:,:), allocatable                :: rand_matrix

    call set_verts_and_coordinate_limits(verts)

    allocate(rand_matrix(5,in%num_particles,in%n_species))
    call RANDOM_NUMBER(rand_matrix)

    call allocate_start_type
    call set_starting_positions(rand_matrix)
    call set_weights
    call set_rest_of_start_type(rand_matrix)

end subroutine calc_starting_conditions

subroutine set_verts_and_coordinate_limits(verts)

    use tetra_physics_mod, only: coord_system, tetra_physics
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, nvert
    use tetra_grid_settings_mod, only: grid_size
    use magdata_in_symfluxcoor_mod, only : raxis,zaxis
    use gorilla_applets_types_mod, only: g

    real(dp), dimension(:,:), allocatable, intent(out)     :: verts
    integer :: i

    g%ind_a = 1 !(R in cylindrical coordinates, s in flux coordinates)
    g%ind_b = 2 !(phi in cylindrical and flux coordinates)
    g%ind_c = 3 !(z in cylindrical coordinates, theta in flux coordinates)
    if (coord_system.eq.2) then
        g%ind_b = 3
        g%ind_c = 2
    endif

    allocate(verts(3,nvert))
    if (coord_system.eq.1) verts = verts_rphiz
    if (coord_system.eq.2) verts = verts_sthetaphi

    g%amin = minval(verts(g%ind_a,:))
    g%amax = maxval(verts(g%ind_a,:))
    g%cmin = minval(verts(g%ind_c,:))
    g%cmax = maxval(verts(g%ind_c,:))

    g%raxis = raxis
    g%zaxis = zaxis

    g%dist_from_o_point_within_grid = 0.0_dp
    do i = 1,3*grid_size(3)
        g%dist_from_o_point_within_grid = max(g%dist_from_o_point_within_grid, &
                                              1.1_dp*sqrt((tetra_physics(i)%x1(1)-raxis)**2 + (tetra_physics(i)%x1(3)-zaxis)**2))
    enddo

end subroutine set_verts_and_coordinate_limits

subroutine allocate_start_type

    use gorilla_applets_types_mod, only: start, in
    use gorilla_applets_settings_mod, only: i_option

    allocate(start%x(3,in%num_particles,in%n_species))
    allocate(start%pitch(in%num_particles,in%n_species))
    allocate(start%energy(in%num_particles,in%n_species))
    allocate(start%weight(in%num_particles,in%n_species))
    allocate(start%jperp(in%num_particles,in%n_species))
    allocate(start%lost(in%num_particles,in%n_species))
    allocate(start%particle_charge(in%n_species))
    allocate(start%particle_mass(in%n_species))
    allocate(start%cm_over_e(in%n_species))
    allocate(start%t(in%n_species))

end subroutine allocate_start_type

subroutine set_starting_positions(rand_matrix)

    use gorilla_applets_types_mod, only: in, start, g
    use tetra_physics_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind
    use constants, only: pi

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix

    !compute starting conditions
    if (in%boole_point_source) then
        if (grid_kind.eq.2) then
            start%x(1,:,:) = 160.0_dp !170.8509699_dp
            start%x(2,:,:) = 0.01_dp
            start%x(3,:,:) = 70.0_dp !8.922304_dp
        elseif (grid_kind.eq.4) then
            start%x(1,:,:) = 205.0_dp
            start%x(2,:,:) = 0.0_dp
            start%x(3,:,:) = 0.0_dp
        endif
        if (coord_system.eq.2) print*, 'error: point source is only implemented for cylindrical coordinate system'
    else
        start%x(g%ind_a,:,:) = g%amin + (g%amax - g%amin)*rand_matrix(g%ind_a,:,:) !r in cylindrical, s in flux coordinates
        start%x(g%ind_b,:,:) = 2*pi*rand_matrix(g%ind_b,:,:) !phi in cylindrical and flux coordinates
        start%x(g%ind_c,:,:) = g%cmin + (g%cmax - g%cmin)*rand_matrix(g%ind_c,:,:) !z in cylindrical, theta in flux coordinates
    endif

end subroutine set_starting_positions

subroutine set_weights

    use gorilla_applets_types_mod, only: in, start, g, c
    use constants, only: pi

    start%weight = in%density*(g%amax-g%amin)*(g%cmax-g%cmin)*2*pi

    if (in%boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start%weight =  start%weight*10/sqrt(pi)
    endif

    c%weight_factor = 1/(start%weight(1,1)*g%amax)

end subroutine set_weights

subroutine set_rest_of_start_type(rand_matrix)

    use gorilla_applets_types_mod, only: in, start
    use tetra_physics_mod, only: cm_over_e, particle_charge, particle_mass
    use gorilla_applets_settings_mod, only: i_option
    use constants, only: echarge,ame,clight
    use constants, only: ev2erg

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix

    start%pitch(:,:) = 2*rand_matrix(4,:,:)-1 !pitch parameter
    start%energy = in%energy_eV
    if (in%boole_boltzmann_energies) then
        start%energy = 5*in%energy_eV*rand_matrix(5,:,:) !boltzmann energy distribution
    endif
    
    if (in%boole_antithetic_variate) then
        start%x(:,1:in%num_particles:2,:) = start%x(:,2:in%num_particles:2,:)
        start%pitch(1:in%num_particles:2,:) = -start%pitch(2:in%num_particles:2,:)
        start%energy(1:in%num_particles:2,:) = start%energy(2:in%num_particles:2,:)
    endif

    start%particle_charge = particle_charge
    start%particle_mass = particle_mass
    start%cm_over_e = cm_over_e
    start%t = in%time_step
    if (i_option.eq.12) then
        start%particle_charge(2) = echarge
        start%particle_mass(2) = ame
        start%cm_over_e(2) = clight*ame/echarge
        start%t(2) = in%time_step/42.0_dp
    endif

    start%v0 = sqrt(2.0_dp*in%energy_eV*ev2erg/start%particle_mass)

end subroutine set_rest_of_start_type

subroutine initialize_exit_data

    use gorilla_applets_types_mod, only: in, exit_data

    allocate(exit_data%lost(in%num_particles))
    allocate(exit_data%t_confined(in%num_particles))
    allocate(exit_data%x(3,in%num_particles))
    allocate(exit_data%vpar(in%num_particles))
    allocate(exit_data%vperp(in%num_particles))
    allocate(exit_data%phi_0_mappings(in%num_particles))
    allocate(exit_data%integration_step(in%num_particles))

    exit_data%lost = 0
    exit_data%t_confined = 0.0_dp
    exit_data%x = 0.0_dp
    exit_data%vpar = 0.0_dp 
    exit_data%vperp = 0.0_dp 
    exit_data%integration_step = 0
    exit_data%phi_0_mappings = 0


end subroutine initialize_exit_data

subroutine calc_poloidal_flux(verts)

    use gorilla_applets_types_mod, only: flux
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: ntetr, tetra_grid
    
    real(dp), dimension(:,:), intent(in) :: verts
    integer :: i
    
    flux%poloidal_max = 0
    flux%poloidal_min = tetra_physics(1)%Aphi1
    do i = 1, ntetr
        flux%poloidal_max = max(flux%poloidal_max,tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
        & (verts([1,2,3],tetra_grid(i)%ind_knot(4))-verts([1,2,3],tetra_grid(i)%ind_knot(1)))))
        flux%poloidal_min = min(flux%poloidal_min,tetra_physics(i)%Aphi1)
    enddo
    
end subroutine calc_poloidal_flux

subroutine calc_collision_coefficients_for_all_tetrahedra(species_in)

    use gorilla_applets_types_mod, only: in, c, start
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_physics_mod, only: particle_mass,particle_charge
    use constants, only: echarge,amp, ame
    use tetra_grid_settings_mod, only: grid_size
    use collis_ions, only: collis_init
    
    integer, intent(in), optional :: species_in
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat
    integer :: Te_unit, Ti_unit, ne_unit
    integer :: i, j
    integer :: species = 1
    real(dp) :: m0, z0

    if (present(species_in)) species = species_in
    
    c%n = 2 !number of background species
    if (.not.allocated(c%dens_mat))   allocate(c%dens_mat(c%n,ntetr))
    if (.not.allocated(c%temp_mat))   allocate(c%temp_mat(c%n,ntetr))
    if (.not.allocated(c%vpar_mat))   allocate(c%vpar_mat(c%n,ntetr))
    if (.not.allocated(c%efcolf_mat)) allocate(c%efcolf_mat(c%n,ntetr))
    if (.not.allocated(c%velrat_mat)) allocate(c%velrat_mat(c%n,ntetr))
    if (.not.allocated(c%enrat_mat))  allocate(c%enrat_mat(c%n,ntetr))
    if (.not.allocated(c%mass))       allocate(c%mass(c%n))
    if (.not.allocated(c%charge_num)) allocate(c%charge_num(c%n))
    if (.not.allocated(c%dens))       allocate(c%dens(c%n))
    if (.not.allocated(c%temp))       allocate(c%temp(c%n))
    if (.not.allocated(efcolf))       allocate(efcolf(c%n))
    if (.not.allocated(velrat))       allocate(velrat(c%n))
    if (.not.allocated(enrat))        allocate(enrat(c%n))
    c%mass = 0
    c%charge_num = 0
    c%mass(1) = 2*amp
    c%mass(c%n) = ame
    !c%mass(2) = 3*amp
    c%charge_num(1) = 1.0d0
    !c%charge_num(2) = 2
    c%charge_num(c%n) = -1.0d0
    c%vpar_mat = 0 !ask Sergei when this will be needed!!!
    m0 = particle_mass
    z0 = particle_charge/echarge
    
    open(newunit = Te_unit, file = 'background/Te_d.dat')
    read(Te_unit,'(e16.9)') (c%temp_mat(2,i),i=1,ntetr/grid_size(2),3)
    close(Te_unit)
    
    open(newunit = Ti_unit, file = 'background/Ti_d.dat')
    read(Ti_unit,'(e16.9)') (c%temp_mat(1,i),i=1,ntetr/grid_size(2),3)
    close(Ti_unit)
    
    open(newunit = ne_unit, file = 'background/ne_d.dat')
    read(ne_unit,'(e16.9)') (c%dens_mat(1,i),i=1,ntetr/grid_size(2),3)
    c%dens_mat(2,:) = c%dens_mat(1,:)
    close(ne_unit)
    
    do i = 1,grid_size(2)-1 !copy data from first phi slice to all other phi slices
        c%temp_mat(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = c%temp_mat(:,1:ntetr/grid_size(2):3)
        c%dens_mat(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = c%dens_mat(:,1:ntetr/grid_size(2):3)
    enddo
    do i = 1,2 !copy data from first tetrahedron of each triangular prism to the two other ones
        c%temp_mat(:,1+i:ntetr:3) = c%temp_mat(:,1:ntetr:3)
        c%dens_mat(:,1+i:ntetr:3) = c%dens_mat(:,1:ntetr:3)
    enddo

    !for now, use constant background profiles
    do i = 1, c%n
        c%dens_mat(i,:) = sum(c%dens_mat(i,:))/ntetr
        c%temp_mat(i,:) = sum(c%temp_mat(i,:))/ntetr
    enddo

    if (.not.in%boole_preserve_energy_and_momentum_during_collisions) then
        do i = 1, ntetr
            do j = 1,c%n
                c%dens(j) = c%dens_mat(j,i)
                c%temp(j) = c%temp_mat(j,i)
            enddo
            call collis_init(m0,z0,c%mass,c%charge_num,c%dens,c%temp,in%energy_eV,start%v0(species), &
                             efcolf,velrat,enrat)
            c%efcolf_mat(:,i) = efcolf
            c%velrat_mat(:,i) = velrat
            c%enrat_mat(:,i) = enrat
        enddo
    endif
    
    if (in%boole_precalc_collisions) then
        if (.not.allocated(c%randcol)) allocate(c%randcol(in%num_particles,c%randcoli,3,in%n_species))
        call RANDOM_NUMBER(c%randcol)
        !3.464102_dp = sqrt(12), this creates a random number with zero average and unit variance
        c%randcol(:,:,1:3:2,:) =  3.464102_dp*(c%randcol(:,:,1:3:2,:)-0.5_dp)
    endif
end subroutine calc_collision_coefficients_for_all_tetrahedra

subroutine perform_background_density_update(i)

    use gorilla_applets_types_mod, only: c, output
    use gorilla_settings_mod, only: boole_time_Hamiltonian
    use tetra_grid_mod, only: ntetr

    integer, intent(in) :: i
    real(dp) :: r=0.99_dp !under-relaxation factor
    real(dp), dimension(ntetr) :: density_from_particle_trajectories

    if (boole_time_Hamiltonian.eqv..false.) then
        print*, "Error, variable 'boole_time_Hamiltonian' must be set to '.true.' for background density update to work"
        stop
    endif

    density_from_particle_trajectories(1:ntetr:3) = output%prism_moments(1,:)
    density_from_particle_trajectories(2:ntetr:3) = output%prism_moments(1,:)
    density_from_particle_trajectories(3:ntetr:3) = output%prism_moments(1,:)

    c%dens_mat(1,:) = r*c%dens_mat(1,:) + (1-r)*density_from_particle_trajectories

    print*, "background density update ", i, " complete"

end subroutine perform_background_density_update

subroutine initialise_electric_potential_type

    use gorilla_applets_types_mod, only: ep, in
    use tetra_grid_mod, only: ntetr, nvert
    use tetra_grid_settings_mod, only: grid_size

    allocate(ep%rho_prism(ntetr/3))
    allocate(ep%rho_flux_tube(grid_size(1)))
    allocate(ep%rho_vert(nvert))
    allocate(ep%phi_elec_from_rho(nvert))
    allocate(ep%average_phi_elec_from_rho(in%n_electric_potential_updates))

    ep%rho_prism = 0
    ep%rho_flux_tube = 0
    ep%rho_vert = 0
    ep%phi_elec_from_rho = 0
    ep%average_phi_elec_from_rho = 0

end subroutine initialise_electric_potential_type

subroutine perform_electric_potential_update(i)

    use gorilla_applets_types_mod, only: c, output, grid_t, in, ep
    use gorilla_settings_mod, only: boole_time_Hamiltonian, coord_system
    use tetra_grid_mod, only: ntetr, nvert
    use tetra_physics_mod, only: make_tetra_physics, phi_elec
    use field_mod, only: ipert

    integer, intent(in) :: i
    real(dp) :: r=0.99_dp !under-relaxation factor

    if (boole_time_Hamiltonian.eqv..false.) then
        print*, "Error, variable 'boole_time_Hamiltonian' must be set to '.true.' for electric potential update to work"
        stop
    endif

    call calc_phi_elec_from_rho

    phi_elec = r*phi_elec + (1-r)*ep%phi_elec_from_rho !make it a module variable there

    !giving i_option = 12 as input variable prevents make_tetra_physics from overwriting phi_elec
    call make_tetra_physics(coord_system,ipert,i_option = 12) 

    ep%average_phi_elec_from_rho(i) = sum(ep%phi_elec_from_rho)/size(ep%phi_elec_from_rho)
    print*, "electric potential update ", i, " complete. Average Delta Phi is ", ep%average_phi_elec_from_rho(i)

end subroutine perform_electric_potential_update

subroutine calc_phi_elec_from_rho

    use gorilla_applets_types_mod, only: in, ep

    real(dp) :: factor = 0.0_dp

    if (in%update_dimension.eq.1) then

        call calc_average_density_per_flux_tube
        call calc_rho_on_vertices
        ep%phi_elec_from_rho = factor * ep%rho_vert
        
    endif

end subroutine calc_phi_elec_from_rho

subroutine calc_average_density_per_flux_tube

    use gorilla_applets_types_mod, only:  g, output, ep
    use tetra_grid_mod, only: nvert
    use tetra_grid_settings_mod, only: grid_size

    integer :: ns, j, ind_prism
    real(dp) :: prism_volumes

    do ns = 1,grid_size(1)
        prism_volumes = 0
        do j = 1, grid_size(2)*grid_size(3)*2
            ind_prism = g%prisms_per_flux_tube(ns,j)
            ep%rho_flux_tube(ns) = ep%rho_flux_tube(ns) + output%prism_moments(1,ind_prism)*output%refined_prism_volumes(ind_prism)
            prism_volumes = prism_volumes + output%refined_prism_volumes(ind_prism)
        enddo
        ep%rho_flux_tube(ns) = ep%rho_flux_tube(ns)/prism_volumes
    enddo

end subroutine calc_average_density_per_flux_tube

subroutine calc_rho_on_vertices

    use tetra_grid_mod, only: nvert
    use tetra_grid_settings_mod, only: grid_size
    use gorilla_applets_types_mod, only: g, ep

    integer :: ns
    real(dp) :: value_to_be_set

    call fill_vector_parts_with_value(ep%rho_vert, g%vertices_per_flux_surface(1,:), ep%rho_flux_tube(1))

    do ns = 2,grid_size(1)
        value_to_be_set = 0.5_dp*(ep%rho_flux_tube(ns-1)+ep%rho_flux_tube(ns))
        call fill_vector_parts_with_value(ep%rho_vert, g%vertices_per_flux_surface(ns,:), value_to_be_set)
    enddo

    call fill_vector_parts_with_value(ep%rho_vert, g%vertices_per_flux_surface(grid_size(1)+1,:), ep%rho_flux_tube(grid_size(1)))

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
        do nphi = 1,grid_size(3)
            do ntheta = 1,grid_size(2)
                ind_vert = (nphi-1)*(grid_size(1)+1)*grid_size(3) + (ns-1)*grid_size(3) + ntheta
                g%vertices_per_flux_surface(ns,i) = ind_vert
                i = i + 1
            enddo
        enddo
    enddo

    !Fill prisms_per_flux_tube
    do ns = 1,grid_size(1)
        i = 1
        do nphi = 1,grid_size(3)
            do ntheta = 1,grid_size(2)
                do j = 1,2 !count through prisms per hexahedron
                    ind_prism = (nphi-1)*grid_size(1)*grid_size(3)*2 + (ns-1)*grid_size(3)*2 + (ntheta-1)*2 + j
                    g%prisms_per_flux_tube(ns,i) = ind_prism
                    i = i + 1
                enddo
            enddo
        enddo
    enddo

end subroutine associate_flux_labels_with_tetrahedra_and_vertices

subroutine prepare_next_round_of_parallelised_particle_pushing(species)

    use gorilla_applets_types_mod, only: start, output, moment_specs
    use constants, only: echarge,ame,clight
    use tetra_physics_mod, only: cm_over_e, particle_charge, particle_mass
    use gorilla_applets_settings_mod, only: i_option

    integer, intent(in), optional :: species

    if (i_option.eq.12) then
        if (.not.present(species)) then
            print*, 'Error: species must be present if i_option = 12. Program terminated.'
            stop
        else
            particle_charge = start%particle_charge(species)
            particle_mass = start%particle_mass(species)
            cm_over_e = start%cm_over_e(species)
        endif
    endif

    call set_weights !weights need to be set again because in orbit_timestep_gorilla_boltzmannn they are multiplied with 
    !some factor, so we need to get rid of it again. Potentially set weights to updated densities
    output%tetr_moments = 0.0_dp
    output%prism_moments = 0.0_dp
    if (moment_specs%boole_squared_moments) output%prism_moments_squared = 0.0_dp

end subroutine prepare_next_round_of_parallelised_particle_pushing

subroutine normalise_prism_moments_and_prism_moments_squared(species)

    use gorilla_applets_types_mod, only: moment_specs, output, in, start
    
    integer, intent(in), optional :: species
    integer :: n
    real(dp) :: time

    time = in%time_step
    if (present(species)) time = start%t(species)
    
    do n = 1,moment_specs%n_moments
        output%prism_moments(n,:) = output%prism_moments(n,:)/(output%prism_volumes*time*in%n_particles)
        if (moment_specs%boole_squared_moments) then
            output%prism_moments_squared(n,:) = output%prism_moments_squared(n,:)/ &
                    (output%prism_volumes**2*time**2*in%n_particles)
        endif
        if (in%boole_refined_sqrt_g) then
            output%prism_moments(n,:) = output%prism_moments(n,:)*output%prism_volumes/output%refined_prism_volumes
            if (moment_specs%boole_squared_moments) then
                output%prism_moments_squared(n,:) = output%prism_moments_squared(n,:)* &
                    output%prism_volumes**2/output%refined_prism_volumes**2
            endif
        endif
    enddo
    
end subroutine normalise_prism_moments_and_prism_moments_squared

subroutine fourier_transform_moments

    use constants, only: pi
    use tetra_grid_settings_mod, only: grid_size
    use gorilla_applets_types_mod, only: moment_specs, output
    use tetra_grid_mod, only: ntetr
    
    integer                                     :: n,m,j,k,p,q,l
    complex                                     :: i
    complex, dimension(:,:,:), allocatable      :: prism_moments_ordered_for_ft
    integer :: n_prisms
    
    n_prisms = ntetr/3
    
    print*, 'fourier transform started'
    
    moment_specs%n_triangles = n_prisms/grid_size(2)
    allocate(prism_moments_ordered_for_ft(moment_specs%n_moments,moment_specs%n_triangles,grid_size(2)))
    do q = 1,grid_size(2)
        prism_moments_ordered_for_ft(:,:,q) = output%prism_moments(:,moment_specs%n_triangles*(q-1)+1:moment_specs%n_triangles*q)
    enddo
    
    if (moment_specs%n_fourier_modes.gt.grid_size(2)) then
        print*, 'moment_specs%n_fourier_modes was chosen to be bigger than n_phi, it is therefore reduced to n_phi'
        moment_specs%n_fourier_modes = grid_size(2)
    endif
    
    i = (0,1)
    
    do p = 1,moment_specs%n_triangles
        do n = 1,moment_specs%n_moments
            do k = 0,moment_specs%n_fourier_modes-1
                do j = 0,grid_size(2)-1
                    output%moments_in_frequency_space(n,p,k+1) = output%moments_in_frequency_space(n,p,k+1) + &
                    1/dble(grid_size(2))*exp(-2*pi*i*j*k/grid_size(2))*prism_moments_ordered_for_ft(n,p,j+1)
                enddo
            enddo
        enddo
    enddo
    
    deallocate(prism_moments_ordered_for_ft)
    
    print*, 'fourier transform finished'
    
end subroutine fourier_transform_moments

subroutine find_minimal_angle_between_curlA_and_tetrahedron_faces
    !tetra_physics(i)%curlA contains contravariant components of the curl of A, so to make them physical components 
    !one has to multiply the second component by R and leave the rest as it is
    !anorm vectors result from the cross product of the vectors in cylindrical coordinates connecting tetrahedron vertices
    !Thus, the second component is physical (containing only r and z components) whereas the first and the third component have to
    !be multiplied by R to make them physical

    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_settings_mod, only: sfc_s_min

    real(dp)                                     :: temp, temp1, R, curlA_norm, normalisation
    real(dp), dimension(3,4)                     :: anorm !normal vectors to tetrahedron faces in physical units
    real(dp), dimension(4)                       :: anormnorm !magnitude of normal vectors
    real(dp), dimension(3)                       :: curlA !physical components of vector potential
    integer                                      :: i, j

    temp = 1.0_dp
    do i = 1,ntetr
        temp1 = 1.0_dp
        anorm = tetra_physics(i)%anorm
        curlA = tetra_physics(i)%curlA !This expression does not contain the metric determinant, but the latter only scales the
        !vector linearly and we are only interested in the normalised vector anyways
        curlA(2) = curlA(2)*verts_rphiz(1,tetra_grid(i)%ind_knot(1))
        curlA_norm = sqrt(curlA(1)**2 + curlA(2)**2 + curlA(3)**2)
        do j = 1,4
            R = verts_rphiz(1,tetra_grid(i)%ind_knot(j))
            anorm(1:3:2,j) = anorm(1:3:2,j)*R
            anormnorm(j) = sqrt(anorm(1,j)**2 + anorm(2,j)**2 + anorm(3,j)**2)
            normalisation = 1/(curlA_norm*anormnorm(j))
            temp1 = min(temp1, abs((anorm(1,j)*curlA(1) + anorm(2,j)*curlA(2) + anorm(3,j)*curlA(3))*normalisation))
        enddo
        ! if (temp1.lt.0.1_dp) then
        !      print*, i, temp1
        !      print*, anormnorm
        !      print*, 'anorm1 = ', anorm(:,1)/anormnorm(1)
        !      print*, 'anorm2 = ', anorm(:,2)/anormnorm(2)
        !      print*, 'anorm3 = ', anorm(:,3)/anormnorm(3)
        !      print*, 'anorm4 = ', anorm(:,4)/anormnorm(4)
        !      print*, 'h = ', tetra_physics(i)%h1_1, tetra_physics(i)%h2_1/tetra_physics(i)%x1(1), tetra_physics(i)%h3_1
        ! endif
        temp = min(temp,temp1)
    enddo
    temp = asin(temp)
    print*, "The minimal angle between the curl of the vector potential and any adjacent tetrahedron face is ", &
            temp, " radiants. sfc_s_min = ", sfc_s_min
end subroutine find_minimal_angle_between_curlA_and_tetrahedron_faces

subroutine analyse_particle_weight_distribution

    use gorilla_applets_types_mod, only: in, start

    integer  :: i
    real(dp) :: maximum_weight, minimum_weight, average_weight

    maximum_weight = start%weight(1,1)
    minimum_weight = start%weight(1,1)
    average_weight = 0

    do i = 1,int(in%n_particles)
        average_weight = average_weight + start%weight(i,1)
        maximum_weight = max(maximum_weight,start%weight(i,1))
        minimum_weight = min(minimum_weight,start%weight(i,1))
    enddo

    average_weight = average_weight/in%n_particles

    print*, "maximum particle weight = ", maximum_weight
    print*, "minimum particle weight = ", minimum_weight
    print*, "average particle weight = ", average_weight

end subroutine analyse_particle_weight_distribution

end module utils_data_pre_and_post_processing_mod