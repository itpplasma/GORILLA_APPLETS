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
    use boltzmann_types_mod, only: moment_specs, in
    
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
    use boltzmann_types_mod, only: moment_specs, output
    
    integer :: n_prisms
    
    n_prisms = ntetr/3
    
    allocate(output%prism_volumes(n_prisms))
    allocate(output%refined_prism_volumes(n_prisms))
    allocate(output%electric_potential(n_prisms))
    allocate(output%boltzmann_density(n_prisms))
    allocate(output%tetr_moments(moment_specs%n_moments,ntetr))
    allocate(output%prism_moments(moment_specs%n_moments,n_prisms))
    if (moment_specs%boole_squared_moments) allocate(output%prism_moments_squared(moment_specs%n_moments,n_prisms))
    allocate(output%moments_in_frequency_space(moment_specs%n_moments,moment_specs%n_triangles,moment_specs%n_fourier_modes))
    
    output%prism_volumes = 0
    output%refined_prism_volumes = 0
    output%electric_potential = 0
    output%boltzmann_density = 0
    output%tetr_moments = 0
    output%prism_moments = 0
    if (moment_specs%boole_squared_moments) output%prism_moments_squared = 0
    output%moments_in_frequency_space = 0
    
end subroutine initialise_output

subroutine calc_starting_conditions(verts)

    use boltzmann_types_mod, only: in
    
    real(dp), dimension(:,:), allocatable, intent(out)     :: verts
    real(dp), dimension(:,:), allocatable                  :: rand_matrix

    call set_verts_and_coordinate_limits(verts)
    allocate(rand_matrix(5,in%num_particles))
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
    use boltzmann_types_mod, only: g

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

    use boltzmann_types_mod, only: start, in

    allocate(start%x(3,in%num_particles))
    allocate(start%pitch(in%num_particles))
    allocate(start%energy(in%num_particles))
    allocate(start%weight(in%num_particles))
    allocate(start%jperp(in%num_particles))

end subroutine allocate_start_type

subroutine set_starting_positions(rand_matrix)

    use boltzmann_types_mod, only: in, start, g
    use tetra_physics_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind
    use constants, only: pi

    real(dp), dimension(:,:), intent(in) :: rand_matrix

    !compute starting conditions
    if (in%boole_point_source) then
        if (grid_kind.eq.2) then
            start%x(1,:) = 160 !170.8509699_dp
            start%x(2,:) = 0.0_dp
            start%x(3,:) = 70 !8.922304_dp
        elseif (grid_kind.eq.4) then
            start%x(1,:) = 205_dp
            start%x(2,:) = 0.0_dp
            start%x(3,:) = 0.0_dp
        endif
        if (coord_system.eq.2) print*, 'error: point source is only implemented for cylindrical coordinate system'
    else
        start%x(g%ind_a,:) = g%amin + (g%amax - g%amin)*rand_matrix(g%ind_a,:) !r in cylindrical, s in flux coordinates
        start%x(g%ind_b,:) = 2*pi*rand_matrix(g%ind_b,:) !phi in cylindrical and flux coordinates
        start%x(g%ind_c,:) = g%cmin + (g%cmax - g%cmin)*rand_matrix(g%ind_c,:) !z in cylindrical, theta in flux coordinates
    endif

end subroutine set_starting_positions

subroutine set_weights

    use boltzmann_types_mod, only: in, start, g, c
    use constants, only: pi

    start%weight = in%density*(g%amax-g%amin)*(g%cmax-g%cmin)*2*pi

    if (in%boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start%weight =  start%weight*10/sqrt(pi)
    endif

    c%weight_factor = 1/(start%weight(1)*g%amax)

end subroutine set_weights

subroutine set_rest_of_start_type(rand_matrix)

    use boltzmann_types_mod, only: in, start

    real(dp), dimension(:,:), intent(in) :: rand_matrix

    start%pitch(:) = 2*rand_matrix(4,:)-1 !pitch parameter
    start%energy = in%energy_eV
    if (in%boole_boltzmann_energies) then
        start%energy = 5*in%energy_eV*rand_matrix(5,:) !boltzmann energy distribution
    endif
    
    if (in%boole_antithetic_variate) then
        start%x(:,1:in%num_particles:2) = start%x(:,2:in%num_particles:2)
        start%pitch(1:in%num_particles:2) = -start%pitch(2:in%num_particles:2)
        start%energy(1:in%num_particles:2) = start%energy(2:in%num_particles:2)
    endif

end subroutine set_rest_of_start_type

subroutine initialize_exit_data

    use boltzmann_types_mod, only: in, exit_data

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

    use boltzmann_types_mod, only: pflux
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: ntetr, tetra_grid
    
    real(dp), dimension(:,:), intent(in) :: verts
    integer :: i
    
    pflux%max = 0
    pflux%min = tetra_physics(1)%Aphi1
    do i = 1, ntetr
        !pflux%max = max(pflux%max,tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
        !& (verts([1,2,3],tetra_grid(i)%ind_knot(4))-verts([1,2,3],tetra_grid(i)%ind_knot(1)))))
        pflux%min = min(pflux%min,tetra_physics(i)%Aphi1)
    enddo
    
end subroutine calc_poloidal_flux

subroutine calc_collision_coefficients_for_all_tetrahedra(v0)

    use boltzmann_types_mod, only: in, c
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_physics_mod, only: particle_mass,particle_charge
    use constants, only: echarge,amp
    use tetra_grid_settings_mod, only: grid_size
    use collis_ions, only: collis_init
    
    real(dp), intent(in) :: v0
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat
    integer :: Te_unit, Ti_unit, ne_unit
    integer :: i, j
    real(dp) :: m0, z0
    
    c%n = 2 !number of background species
    allocate(c%dens_mat(c%n-1,ntetr))
    allocate(c%temp_mat(c%n,ntetr))
    allocate(c%vpar_mat(c%n,ntetr))
    allocate(c%efcolf_mat(c%n,ntetr))
    allocate(c%velrat_mat(c%n,ntetr))
    allocate(c%enrat_mat(c%n,ntetr))
    allocate(c%mass_num(c%n-1))
    allocate(c%charge_num(c%n-1))
    allocate(c%dens(c%n))
    allocate(c%temp(c%n))
    allocate(efcolf(c%n))
    allocate(velrat(c%n))
    allocate(enrat(c%n))
    c%mass_num = 0
    c%charge_num = 0
    c%mass_num(1) = 2
    !c%mass_num(2) = 3
    c%charge_num(1) = 1
    !c%charge_num(2) = 2
    c%vpar_mat = 0 !ask Sergei when this will be needed!!!
    m0 = particle_mass/amp
    z0 = particle_charge/echarge
    
    open(newunit = Te_unit, file = 'background/Te_d.dat')
    read(Te_unit,'(e16.9)') (c%temp_mat(2,i),i=1,ntetr/grid_size(2),3)
    close(Te_unit)
    
    open(newunit = Ti_unit, file = 'background/Ti_d.dat')
    read(Ti_unit,'(e16.9)') (c%temp_mat(1,i),i=1,ntetr/grid_size(2),3)
    close(Ti_unit)
    
    open(newunit = ne_unit, file = 'background/ne_d.dat')
    read(ne_unit,'(e16.9)') (c%dens_mat(1,i),i=1,ntetr/grid_size(2),3)
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
    do i = 1, c%n-1
        c%dens_mat(i,:) = sum(c%dens_mat(i,:))/ntetr
        c%temp_mat(i,:) = sum(c%temp_mat(i,:))/ntetr
    enddo
    c%temp_mat(c%n,:) = sum(c%temp_mat(c%n,:))/ntetr

    if (.not.in%boole_preserve_energy_and_momentum_during_collisions) then
        do i = 1, ntetr
            do j = 1,c%n
                !> if statement because electron density will be calculated in collis init
                if (j.lt.c%n) c%dens(j) = c%dens_mat(j,i)
                c%temp(j) = c%temp_mat(j,i)
            enddo
            call collis_init(m0,z0,c%mass_num,c%charge_num,c%dens,c%temp,in%energy_eV,v0,efcolf,velrat,enrat)
            c%efcolf_mat(:,i) = efcolf
            c%velrat_mat(:,i) = velrat
            c%enrat_mat(:,i) = enrat
        enddo
    endif
    
    if (in%boole_precalc_collisions) then
        allocate(c%randcol(in%num_particles,c%randcoli,3))
        call RANDOM_NUMBER(c%randcol)
        !3.464102_dp = sqrt(12), this creates a random number with zero average and unit variance
        c%randcol(:,:,1:3:2) =  3.464102_dp*(c%randcol(:,:,1:3:2)-0.5_dp)
    endif
end subroutine calc_collision_coefficients_for_all_tetrahedra

subroutine perform_background_density_update(i)

    use boltzmann_types_mod, only: c, output
    use gorilla_settings_mod, only: boole_time_Hamiltonian

    integer, intent(in) :: i
    real(dp) :: r=0.99_dp

    if (boole_time_Hamiltonian.eqv..false.) then
        print*, "Error, variable 'boole_time_Hamiltonian' must be set to '.true.' for background density update to work"
        stop
    endif

    c%dens_mat(1,:) = r*c%dens_mat(1,:) + (1-r)*output%tetr_moments(1,:)

    print*, "background density update ", i, " complete"

end subroutine perform_background_density_update

subroutine normalise_prism_moments_and_prism_moments_squared

    use boltzmann_types_mod, only: moment_specs, output, in
    
    integer :: n
    
    do n = 1,moment_specs%n_moments
        output%prism_moments(n,:) = output%prism_moments(n,:)/(output%prism_volumes*in%time_step*in%n_particles)
        if (moment_specs%boole_squared_moments) then
            output%prism_moments_squared(n,:) = output%prism_moments_squared(n,:)/ &
                    (output%prism_volumes**2*in%time_step**2*in%n_particles)
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
    use boltzmann_types_mod, only: moment_specs, output
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

    use boltzmann_types_mod, only: in, start

    integer  :: i
    real(dp) :: maximum_weight, minimum_weight, average_weight

    maximum_weight = start%weight(1)
    minimum_weight = start%weight(1)
    average_weight = 0

    do i = 1,int(in%n_particles)
        average_weight = average_weight + start%weight(i)
        maximum_weight = max(maximum_weight,start%weight(i))
        minimum_weight = min(minimum_weight,start%weight(i))
    enddo

    average_weight = average_weight/in%n_particles

    print*, "maximum particle weight = ", maximum_weight
    print*, "minimum particle weight = ", minimum_weight
    print*, "average particle weight = ", average_weight

end subroutine analyse_particle_weight_distribution

end module utils_data_pre_and_post_processing_mod