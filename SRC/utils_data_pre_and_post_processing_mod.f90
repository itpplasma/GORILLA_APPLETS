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
    use gorilla_applets_types_mod, only: moment_specs, output, in
    use tetra_grid_settings_mod, only: grid_size
    
    integer :: n_prisms
    
    n_prisms = ntetr/3
    
    allocate(output%prism_volumes(n_prisms))
    allocate(output%refined_prism_volumes(n_prisms))
    allocate(output%electric_potential(n_prisms))
    allocate(output%boltzmann_density(n_prisms))
    allocate(output%radial_flux(grid_size(1)+1))
    allocate(output%tetr_moments(moment_specs%n_moments,ntetr,in%n_species))
    allocate(output%prism_moments(moment_specs%n_moments,n_prisms,in%n_species))
    if (moment_specs%boole_squared_moments) allocate(output%prism_moments_squared(moment_specs%n_moments,n_prisms,in%n_species))
    allocate(output%moments_in_frequency_space(moment_specs%n_moments,moment_specs%n_triangles,moment_specs%n_fourier_modes))
    
    output%prism_volumes = 0.0_dp
    output%refined_prism_volumes = 0.0_dp
    output%electric_potential = 0.0_dp
    output%boltzmann_density = 0.0_dp
    output%radial_flux = 0.0_dp
    output%tetr_moments = 0.0_dp
    output%prism_moments = 0.0_dp
    if (moment_specs%boole_squared_moments) output%prism_moments_squared = 0.0_dp
    output%moments_in_frequency_space = 0.0_dp
    
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
        start%particle_charge(2) = -echarge
        start%particle_mass(2) = ame
        start%cm_over_e(2) = -clight*ame/echarge
        start%t(2) = in%time_step/42.0_dp
    endif

    start%v0 = sqrt(2.0_dp*in%energy_eV*ev2erg/start%particle_mass)

end subroutine set_rest_of_start_type

subroutine initialize_exit_data(n_particles_in)

    use gorilla_applets_types_mod, only: in, exit_data

    integer, intent(in), optional :: n_particles_in 
    integer                       :: n_particles

    if (allocated(exit_data%lost))              deallocate(exit_data%lost)
    if (allocated(exit_data%t_confined))        deallocate(exit_data%t_confined)
    if (allocated(exit_data%x))                 deallocate(exit_data%x)
    if (allocated(exit_data%vpar))              deallocate(exit_data%vpar)
    if (allocated(exit_data%vperp))             deallocate(exit_data%vperp)
    if (allocated(exit_data%phi_0_mappings))    deallocate(exit_data%phi_0_mappings)
    if (allocated(exit_data%integration_step))  deallocate(exit_data%integration_step)

    if(present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

    allocate(exit_data%lost(n_particles,in%n_species))
    allocate(exit_data%t_confined(n_particles,in%n_species))
    allocate(exit_data%x(3,n_particles,in%n_species))
    allocate(exit_data%vpar(n_particles,in%n_species))
    allocate(exit_data%vperp(n_particles,in%n_species))
    allocate(exit_data%phi_0_mappings(n_particles,in%n_species))
    allocate(exit_data%integration_step(n_particles,in%n_species))

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

    use gorilla_applets_types_mod, only: in, c, start, s
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_physics_mod, only: particle_mass,particle_charge, tetra_physics
    use constants, only: echarge,amp, ame, ev2erg
    use tetra_grid_settings_mod, only: grid_size
    use collis_ions, only: collis_init
    use gorilla_settings_mod, only: coord_system
    
    integer, intent(in), optional :: species_in
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat
    integer :: Te_unit, Ti_unit, ne_unit
    integer :: i, j
    integer :: species = 1
    real(dp) :: m0, z0, n0, s_value, v0

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

    !If working in flux coordinates, use background density linear in s
    if (coord_system.eq.2) then
        n0 = 3.0_dp * 10.0_dp**13
        do i = 1,ntetr/grid_size(2),3
            s_value = tetra_physics(i)%x1(1)
            c%dens_mat(1,i) = n0*(1-s_value*0.9_dp)
        enddo
        c%dens_mat(2,:) = c%dens_mat(1,:)
    endif
    
    do i = 1,grid_size(2)-1 !copy data from first phi slice to all other phi slices
        c%temp_mat(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = c%temp_mat(:,1:ntetr/grid_size(2):3)
        c%dens_mat(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = c%dens_mat(:,1:ntetr/grid_size(2):3)
    enddo
    do i = 1,2 !copy data from first tetrahedron of each triangular prism to the two other ones
        c%temp_mat(:,1+i:ntetr:3) = c%temp_mat(:,1:ntetr:3)
        c%dens_mat(:,1+i:ntetr:3) = c%dens_mat(:,1:ntetr:3)
    enddo

    !for now, use constant background temperature
    do i = 1, c%n
        c%temp_mat(i,:) = sum(c%temp_mat(i,:))/ntetr
        if (coord_system.eq.1) c%dens_mat(i,:) = sum(c%dens_mat(i,:))/ntetr
    enddo
    
    if (coord_system.eq.2) c%temp_mat(2,:) = s%temperature
    if (coord_system.eq.2) c%temp_mat(1,:) = in%energy_eV
    if (coord_system.eq.2) c%dens_mat = in%density


    if (.not.in%boole_preserve_energy_and_momentum_during_collisions) then
        do i = 1, ntetr
            do j = 1,c%n
                c%dens(j) = c%dens_mat(j,i)
                c%temp(j) = c%temp_mat(j,i)
            enddo
            call collis_init(m0,z0,c%mass,c%charge_num,c%dens,c%temp,in%energy_eV,v0, &
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

    density_from_particle_trajectories(1:ntetr:3) = real(output%prism_moments(1,:,1))
    density_from_particle_trajectories(2:ntetr:3) = real(output%prism_moments(1,:,1))
    density_from_particle_trajectories(3:ntetr:3) = real(output%prism_moments(1,:,1))

    c%dens_mat(1,:) = r*c%dens_mat(1,:) + (1-r)*density_from_particle_trajectories

    print*, "background density update ", i, " complete"

end subroutine perform_background_density_update

subroutine prepare_next_round_of_parallelised_particle_pushing(species)

    use gorilla_applets_types_mod, only: start, output, moment_specs, counter
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

    output%tetr_moments(:,:,species) = 0.0_dp
    output%prism_moments(:,:,species) = 0.0_dp
    if (moment_specs%boole_squared_moments) output%prism_moments_squared(:,:,species) = 0.0_dp

    counter%lost_particles = 0
    counter%lost_inside = 0
    counter%tetr_pushings = 0
    counter%phi_0_mappings = 0
    counter%integration_steps = 0

end subroutine prepare_next_round_of_parallelised_particle_pushing

subroutine normalise_prism_moments_and_prism_moments_squared(species_in)

    use gorilla_applets_types_mod, only: moment_specs, output, in, start
    
    integer, intent(in), optional :: species_in
    integer :: species = 1
    integer :: n
    real(dp) :: time

    time = in%time_step
    if (present(species_in)) species = species_in
    time = start%t(species)
    
    do n = 1,moment_specs%n_moments
        output%prism_moments(n,:,species) = output%prism_moments(n,:,species)/(output%prism_volumes*time*in%n_particles)
        if (moment_specs%boole_squared_moments) then
            output%prism_moments_squared(n,:,species) = output%prism_moments_squared(n,:,species)/ &
                    (output%prism_volumes**2*time**2*in%n_particles)
        endif
        if (in%boole_refined_sqrt_g) then
            output%prism_moments(n,:,species) = output%prism_moments(n,:,species)*output%prism_volumes/output%refined_prism_volumes
            if (moment_specs%boole_squared_moments) then
                output%prism_moments_squared(n,:,species) = output%prism_moments_squared(n,:,species)* &
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
        prism_moments_ordered_for_ft(:,:,q) = output%prism_moments(:,moment_specs%n_triangles*(q-1)+1:moment_specs%n_triangles*q,1)
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