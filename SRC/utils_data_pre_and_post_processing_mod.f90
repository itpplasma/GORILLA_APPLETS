module utils_data_pre_and_post_processing_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64

    implicit none

    integer, dimension(:), allocatable, protected :: master_seed
    integer(int64),                     protected :: master_seed_u64 = 0_int64
    logical,                            protected :: master_seed_loaded = .false.

contains

subroutine set_seed_for_random_numbers

    integer,dimension(:), allocatable   :: seed
    integer                             :: seed_inp_unit
    integer                             :: n, i
    integer(int64)                      :: acc, tmp

    open(newunit = seed_inp_unit, file='seed.inp', status='old',action = 'read')
    read(seed_inp_unit,*) n
    allocate(seed(n))
    read(seed_inp_unit,*) seed
    close(seed_inp_unit)
    CALL RANDOM_SEED (PUT=seed)

    ! Retain the loaded seed and derive a 64-bit master key from it, so that
    ! per-thread seeding inside parallel loops can be a deterministic function
    ! of seed.inp (see initialise_seed_for_random_numbers_for_each_thread).
    if (allocated(master_seed)) deallocate(master_seed)
    allocate(master_seed(n))
    master_seed = seed

    acc = 0_int64
    do i = 1, n
        acc = ieor(acc, int(seed(i), int64))
        call splitmix64_next(acc, tmp)
    enddo
    master_seed_u64 = acc
    master_seed_loaded = .true.

    deallocate(seed)

end subroutine set_seed_for_random_numbers

subroutine splitmix64_next(state, z)
    ! One step of Vigna's SplitMix64 bit mixer. Advances `state` in place and
    ! returns the next well-diffused 64-bit word `z`. Used both to fold the
    ! master seed into a single key and to fan a key out into a per-thread
    ! Fortran RNG state array.

    integer(int64), intent(inout) :: state
    integer(int64), intent(out)   :: z
    integer(int64), parameter     :: gamma = -7046029254386353131_int64  ! 0x9E3779B97F4A7C15
    integer(int64), parameter     :: m1    = -4658895280553007687_int64  ! 0xBF58476D1CE4E5B9
    integer(int64), parameter     :: m2    = -7723592293110705173_int64  ! 0x94D049BB133111EB

    state = state + gamma
    z = state
    z = ieor(z, ishft(z, -30)) * m1
    z = ieor(z, ishft(z, -27)) * m2
    z = ieor(z, ishft(z, -31))

end subroutine splitmix64_next

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

    !Free any pre-existing output arrays so a fresh build can proceed without external bookkeeping.
    call deallocate_output

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

subroutine setup_applet

    !Initialise GORILLA and the moment output, after the mode has read its input file

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert

    call get_ipert()
    call initialize_gorilla(i_option,ipert)

    call set_moment_specifications
    call initialise_output

end subroutine setup_applet

subroutine print_applet_summary(particle_name, particle_mass, v0, particle_charge)

    !Print the counters and particle parameters at the end of a run
    !particle_name ... 'particle' or 'ion'

    use constants, only: ev2erg
    use tetra_grid_settings_mod, only: grid_kind
    use gorilla_applets_types_mod, only: counter, c, in

    character(len=*), intent(in) :: particle_name
    real(dp), intent(in)         :: particle_mass, v0, particle_charge

    if (in%boole_precalc_collisions) print*, "maxcol = ", c%maxcol
    print*, 'Number of lost '//particle_name//'s', counter%lost_particles
    print*, 'average number of pushings = ', counter%tetr_pushings/in%n_particles
    print*, 'average number of toroidal revolutions = ', counter%phi_0_mappings/in%n_particles
    print*, 'average number of integration steps = ', counter%integration_steps/in%n_particles
    print*, 'number of times the pusher gave up on a particle = ', counter%pusher_failures
    print*, particle_name//' mass = ', particle_mass
    print*, 'absolute value of '//particle_name//' velocity = ', v0
    print*, particle_name//' charge = ', particle_charge
    print*, 'temperature = ', ev2erg*in%energy_eV
    print*, 'energy in eV = ', in%energy_eV
    print*, 'tracing time in seconds = ', in%time_step
    if((grid_kind.eq.2).or.(grid_kind.eq.3)) then
         print*, 'number of times that particles were pushed across the inside hole = ', counter%lost_inside
    endif

end subroutine print_applet_summary

subroutine deallocate_output
!
! Deallocates arrays in the output type to allow reallocation when grid parameters change.
!
    use gorilla_applets_types_mod, only: output

    if (allocated(output%prism_volumes)) deallocate(output%prism_volumes)
    if (allocated(output%refined_prism_volumes)) deallocate(output%refined_prism_volumes)
    if (allocated(output%electric_potential)) deallocate(output%electric_potential)
    if (allocated(output%boltzmann_density)) deallocate(output%boltzmann_density)
    if (allocated(output%radial_flux)) deallocate(output%radial_flux)
    if (allocated(output%tetr_moments)) deallocate(output%tetr_moments)
    if (allocated(output%prism_moments)) deallocate(output%prism_moments)
    if (allocated(output%prism_moments_squared)) deallocate(output%prism_moments_squared)
    if (allocated(output%moments_in_frequency_space)) deallocate(output%moments_in_frequency_space)

end subroutine deallocate_output

subroutine calc_starting_conditions(verts)

    use gorilla_applets_types_mod, only: in, start

    real(dp), dimension(:,:), allocatable, intent(out), optional :: verts

    call set_verts_and_coordinate_limits(verts)

    call allocate_start_type
    call allocate_weights
    call set_starting_positions()
    call set_rest_of_start_type()

end subroutine calc_starting_conditions

subroutine set_verts_and_coordinate_limits(verts)

    use tetra_physics_mod, only: coord_system, tetra_physics
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, nvert
    use tetra_grid_settings_mod, only: grid_size
    use magdata_in_symfluxcoor_mod, only : raxis,zaxis
    use gorilla_applets_types_mod, only: g

    real(dp), dimension(:,:), allocatable, intent(out), optional :: verts
    real(dp), dimension(:,:), allocatable :: verts_local
    integer :: i

    g%ind_a = 1 !(R in cylindrical coordinates, s in flux coordinates)
    g%ind_b = 2 !(phi in cylindrical and flux coordinates)
    g%ind_c = 3 !(z in cylindrical coordinates, theta in flux coordinates)
    if (coord_system.eq.2) then
        g%ind_b = 3
        g%ind_c = 2
    endif

    allocate(verts_local(3,nvert))
    if (coord_system.eq.1) verts_local = verts_rphiz
    if (coord_system.eq.2) verts_local = verts_sthetaphi

    g%amin = minval(verts_local(g%ind_a,:))
    g%amax = maxval(verts_local(g%ind_a,:))
    g%cmin = minval(verts_local(g%ind_c,:))
    g%cmax = maxval(verts_local(g%ind_c,:))

    g%raxis = raxis
    g%zaxis = zaxis

    g%dist_from_o_point_within_grid = 0.0_dp
    do i = 1,3*grid_size(3)
        g%dist_from_o_point_within_grid = max(g%dist_from_o_point_within_grid, &
                                              1.1_dp*sqrt((tetra_physics(i)%x1(1)-raxis)**2 + (tetra_physics(i)%x1(3)-zaxis)**2))
    enddo

    if (present(verts)) call move_alloc(verts_local, verts)

end subroutine set_verts_and_coordinate_limits

subroutine allocate_start_type

    use gorilla_applets_types_mod, only: start, in
    use gorilla_applets_settings_mod, only: i_option

    allocate(start%x(3,in%num_particles,in%n_species))
    allocate(start%pitch(in%num_particles,in%n_species))
    allocate(start%energy(in%num_particles,in%n_species))
    allocate(start%jperp(in%num_particles,in%n_species))
    allocate(start%density(in%num_particles,in%n_species))
    allocate(start%lost(in%num_particles,in%n_species))
    allocate(start%particle_charge(in%n_species))
    allocate(start%particle_mass(in%n_species))
    allocate(start%cm_over_e(in%n_species))
    allocate(start%t(in%n_species))
    allocate(start%v0(in%n_species))

    start%density = in%density

end subroutine allocate_start_type

subroutine allocate_weights

    use gorilla_applets_types_mod, only: in, weights

    allocate(weights%w(in%num_particles,in%n_species))
    weights%w = (0.0_dp, 0.0_dp)
    if (in%boole_delta_f) then
        allocate(weights%original(in%num_particles,in%n_species))
        weights%original = (0.0_dp, 0.0_dp)
    end if

end subroutine allocate_weights

subroutine set_starting_positions()

    use gorilla_applets_types_mod, only: in, start, g
    use tetra_physics_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind
    use constants, only: pi
    use marker_distribution_mod, only: pdf_flat, init_distribution_3d, sample_array_3d

    real(dp) :: xmin(3), xmax(3)
    integer :: i_species

    !compute starting conditions
    if (in%boole_point_source) then
        if (grid_kind.eq.2) then
            xmin = [209.0_dp, 0.01_dp, 10.0_dp]
            xmax = xmin
        elseif (grid_kind.eq.4) then
            xmin = [205.0_dp, 0.0_dp, 0.0_dp]
            xmax = xmin
        endif
        if (coord_system.eq.2) print*, 'error: point source is only implemented for cylindrical coordinate system'
    else
        ! Set bounds for uniform sampling
        xmin(g%ind_a) = g%amin
        xmax(g%ind_a) = g%amax
        xmin(g%ind_b) = 0.0_dp
        xmax(g%ind_b) = 2*pi
        xmin(g%ind_c) = g%cmin
        xmax(g%ind_c) = g%cmax
    endif

    ! Initialize position distribution if not already done
    if (.not. start%dist_position%initialized) then
        call init_distribution_3d(start%dist_position, pdf_flat, xmin, xmax)
    endif

    ! Sample positions for all species
    do i_species = 1, in%n_species
        call sample_array_3d(start%dist_position, start%x(:,:,i_species))
    enddo

end subroutine set_starting_positions

subroutine init_velocity_space_distributions()

    !Marker distributions of pitch (uniform on [-1, 1]) and energy (fixed for monoenergetic runs,
    !uniform on [0, 5*T] otherwise), used for sampling and for evaluating the marker weights

    use gorilla_applets_types_mod, only: in, start
    use marker_distribution_mod, only: pdf_flat, init_distribution_1d

    if (.not. start%dist_lambda%initialized) then
        call init_distribution_1d(start%dist_lambda, pdf_flat, -1.0_dp, 1.0_dp)
    endif

    if (.not. start%dist_energy%initialized) then
        if (in%boole_monoenergetic) then
            call init_distribution_1d(start%dist_energy, pdf_flat, in%energy_eV, in%energy_eV)
        else
            call init_distribution_1d(start%dist_energy, pdf_flat, 0.0_dp, 5.0_dp*in%energy_eV)
        endif
    endif

end subroutine init_velocity_space_distributions

subroutine set_rest_of_start_type()

    use gorilla_applets_types_mod, only: in, start
    use marker_distribution_mod, only: sample_array_1d

    integer :: i_species

    call init_velocity_space_distributions

    ! Sample pitch and energy for all species
    do i_species = 1, in%n_species
        call sample_array_1d(start%dist_lambda, in%num_particles, start%pitch(:,i_species))
        call sample_array_1d(start%dist_energy, in%num_particles, start%energy(:,i_species))
    enddo

    if (in%boole_antithetic_variate) then
        start%x(:,1:in%num_particles:2,:) = start%x(:,2:in%num_particles:2,:)
        start%pitch(1:in%num_particles:2,:) = -start%pitch(2:in%num_particles:2,:)
        start%energy(1:in%num_particles:2,:) = start%energy(2:in%num_particles:2,:)
    endif

    call set_species_properties

end subroutine set_rest_of_start_type

subroutine set_species_properties()

    !Per-species charge, mass, tracing time and reference velocity, read by the shared particle
    !pushing and weight routines. Applets placing their markers themselves (divertor) call this
    !directly instead of set_rest_of_start_type.

    use gorilla_applets_types_mod, only: in, start
    use tetra_physics_mod, only: cm_over_e, particle_charge, particle_mass
    use gorilla_applets_settings_mod, only: i_option
    use constants, only: echarge,ame,clight
    use constants, only: ev2erg

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
    start%lost = .false.

end subroutine set_species_properties

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
    if (allocated(exit_data%flux_surface))      deallocate(exit_data%flux_surface)

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
    allocate(exit_data%flux_surface(n_particles,in%n_species))

    exit_data%lost = 0
    exit_data%t_confined = 0.0_dp
    exit_data%x = 0.0_dp
    exit_data%vpar = 0.0_dp
    exit_data%vperp = 0.0_dp
    exit_data%integration_step = 0
    exit_data%phi_0_mappings = 0
    exit_data%flux_surface = 0.0_dp

end subroutine initialize_exit_data

subroutine calc_poloidal_flux(verts)

    use gorilla_applets_types_mod, only: flux
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: ntetr, tetra_grid
    
    real(dp), dimension(:,:), intent(in) :: verts
    integer :: i, k
    real(dp) :: psi_vertex

    !Range of A_phi over all vertices of the grid
    flux%poloidal_min = huge(1.0_dp)
    flux%poloidal_max = -huge(1.0_dp)
    do i = 1, ntetr
        do k = 1, 4
            psi_vertex = tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
            & (verts([1,2,3],tetra_grid(i)%ind_knot(k))-verts([1,2,3],tetra_grid(i)%ind_knot(1))))
            flux%poloidal_min = min(flux%poloidal_min, psi_vertex)
            flux%poloidal_max = max(flux%poloidal_max, psi_vertex)
        enddo
    enddo
    
end subroutine calc_poloidal_flux

elemental function normalised_poloidal_flux(psi) result(s)

    !Poloidal flux psi = A_phi normalised to [0, 1] by its range over the computational domain
    !(flux%poloidal_min/max from calc_poloidal_flux); the radial coordinate in cylindrical coordinates
    !for axisymmetric devices. A_phi has no fixed sign or offset, so dividing by poloidal_max alone
    !can give s < 0.

    use gorilla_applets_types_mod, only: flux

    real(dp), intent(in) :: psi
    real(dp) :: s

    s = (psi - flux%poloidal_min)/(flux%poloidal_max - flux%poloidal_min)

end function normalised_poloidal_flux

subroutine calc_collision_coefficients_for_all_tetrahedra(species_in)

    use gorilla_applets_types_mod, only: in, c, start, s, g
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_physics_mod, only: particle_mass,particle_charge, tetra_physics
    use constants, only: echarge,amp, ame, ev2erg, pi
    use tetra_grid_settings_mod, only: grid_size, grid_kind
    use collis_ions, only: collis_init
    use gorilla_settings_mod, only: coord_system
    use gorilla_applets_settings_mod, only: i_option
    use marker_distribution_mod, only: linear_profile_factor

    integer, intent(in), optional :: species_in
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat
    real(dp), dimension(:), allocatable :: efcolf_l, velrat_l, enrat_l, dens_l, temp_l
    real(dp) :: v0_l
    integer :: i, j
    integer :: species
    real(dp) :: m0, z0, n0, s_value, v0
    logical :: boole_T_and_n_from_files = .false.
    logical :: boole_analytic_temperature_profile

    species = 1

    if (present(species_in)) species = species_in
    ! The analytic SCEF temperature model must not replace other applet profiles.
    boole_analytic_temperature_profile = (i_option == 12)

    call set_c

    !The energy/momentum-conserving collisions update a background reservoir per tetrahedron, whose
    !particle number is density times the physical tetrahedron volume
    if (in%boole_preserve_energy_and_momentum_during_collisions) call calc_physical_tetrahedron_volumes

    if (.not.allocated(efcolf))       allocate(efcolf(c%n))
    if (.not.allocated(velrat))       allocate(velrat(c%n))
    if (.not.allocated(enrat))        allocate(enrat(c%n))

    m0 = particle_mass
    z0 = particle_charge/echarge


    c%temp_mat = in%energy_eV
    if (i_option.eq.12) c%temp_mat(c%n, :) = s%temperature
    ! Ion densities only; electron density is set below via quasi-neutrality.
    c%dens_mat(1:c%n-1, :) = in%density

    if (boole_T_and_n_from_files) call get_T_and_n_from_files

    ! Apply linear density and/or temperature profiles based on input settings
    if (in%boole_linear_density_simulation .or. in%boole_linear_temperature_simulation .or. boole_analytic_temperature_profile) then
        do i = 1, ntetr/grid_size(2), 3
            ! Compute s_value (normalized flux coordinate) for this tetrahedron
            if (coord_system == 2) then
                ! Flux coordinates: use s-coordinate directly
                s_value = tetra_physics(i)%x1(1)
            else if (grid_kind /= 3) then
                ! Cylindrical coordinates with axisymmetric device: use poloidal flux from A_phi
                s_value = normalised_poloidal_flux(tetra_physics(i)%Aphi1)
            else
                ! grid_kind == 3 (stellarator) with cylindrical coordinates: not supported
                print*, 'Error in calc_collision_coefficients_for_all_tetrahedra: Linear profiles with &
                        &cylindrical coordinates only supported for axisymmetric devices.'
                stop
            endif

            if (in%boole_linear_density_simulation) then
                c%dens_mat(:,i) = in%density * linear_profile_factor(s_value)
            endif

            if (in%boole_linear_temperature_simulation) then
                c%temp_mat(:,i) = in%energy_eV * linear_profile_factor(s_value)
            endif

            if (boole_analytic_temperature_profile) then
                c%temp_mat(1,   i) = 1770.0_dp*(0.06_dp + 0.94_dp*(1.0_dp - s_value**2.9_dp)**2)
                c%temp_mat(c%n, i) = c%temp_mat(1, i) + 3190.0_dp*(1.0_dp - s_value**2.1_dp)**5.4_dp
            endif
        enddo
    endif

    do i = 1,grid_size(2)-1 !copy data from first phi slice to all other phi slices
        c%temp_mat(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = c%temp_mat(:,1:ntetr/grid_size(2):3)
        c%dens_mat(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = c%dens_mat(:,1:ntetr/grid_size(2):3)
    enddo
    do i = 1,2 !copy data from first tetrahedron of each triangular prism to the two other ones
        c%temp_mat(:,1+i:ntetr:3) = c%temp_mat(:,1:ntetr:3)
        c%dens_mat(:,1+i:ntetr:3) = c%dens_mat(:,1:ntetr:3)
    enddo

    ! Enforce quasi-neutrality: electron density = sum of Z_i * n_i over ion species.
    c%dens_mat(c%n, :) = matmul(c%charge_num(1:c%n-1), c%dens_mat(1:c%n-1, :))

    if (.not.in%boole_preserve_energy_and_momentum_during_collisions) then
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,efcolf_l,velrat_l,enrat_l,dens_l,temp_l,v0_l)
        allocate(efcolf_l(c%n), velrat_l(c%n), enrat_l(c%n), dens_l(c%n), temp_l(c%n))
        !$OMP DO
        do i = 1, ntetr
            do j = 1,c%n
                dens_l(j) = c%dens_mat(j,i)
                temp_l(j) = c%temp_mat(j,i)
            enddo
            call collis_init(m0,z0,c%mass,c%charge_num,dens_l,temp_l,in%energy_eV,v0_l, &
                             efcolf_l,velrat_l,enrat_l)
            c%efcolf_mat(:,i) = efcolf_l
            c%velrat_mat(:,i) = velrat_l
            c%enrat_mat(:,i) = enrat_l
        enddo
        !$OMP END DO
        deallocate(efcolf_l, velrat_l, enrat_l, dens_l, temp_l)
        !$OMP END PARALLEL
    endif
    
    if (in%boole_precalc_collisions) then
        if (.not.allocated(c%randcol)) allocate(c%randcol(in%num_particles,c%randcoli,3,in%n_species))
        call RANDOM_NUMBER(c%randcol)
        !3.464102_dp = sqrt(12), this creates a random number with zero average and unit variance
        c%randcol(:,:,1:3:2,:) =  3.464102_dp*(c%randcol(:,:,1:3:2,:)-0.5_dp)
    endif

end subroutine calc_collision_coefficients_for_all_tetrahedra

subroutine calc_physical_tetrahedron_volumes

    !c%tetr_volume = coordinate volume times the Jacobian at the centroid, which is exact for the
    !Jacobian R of cylindrical coordinates (linear in the tetrahedron) and uses the linearised
    !sqrt(g) of GORILLA in flux coordinates. Toroidal angles are unwrapped relative to vertex 1.

    use gorilla_applets_types_mod, only: c
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_rphiz, verts_sthetaphi
    use tetra_physics_mod, only: tetra_physics, coord_system
    use constants, only: pi

    integer :: i, k, ind_phi
    real(dp), dimension(3,4) :: v
    real(dp), dimension(3) :: centroid
    real(dp), dimension(3,3) :: edges
    real(dp) :: coordinate_volume, jacobian

    if (allocated(c%tetr_volume)) deallocate(c%tetr_volume)
    allocate(c%tetr_volume(ntetr))

    ind_phi = 2
    if (coord_system.eq.2) ind_phi = 3

    do i = 1, ntetr
        do k = 1, 4
            if (coord_system.eq.1) then
                v(:,k) = verts_rphiz(:,tetra_grid(i)%ind_knot(k))
            else
                v(:,k) = verts_sthetaphi(:,tetra_grid(i)%ind_knot(k))
            endif
        enddo
        do k = 2, 4
            v(ind_phi,k) = v(ind_phi,1) + modulo(v(ind_phi,k) - v(ind_phi,1) + pi, 2.0_dp*pi) - pi
        enddo
        edges = v(:,2:4) - spread(v(:,1), 2, 3)
        coordinate_volume = abs(edges(1,1)*(edges(2,2)*edges(3,3) - edges(3,2)*edges(2,3)) &
                              - edges(1,2)*(edges(2,1)*edges(3,3) - edges(3,1)*edges(2,3)) &
                              + edges(1,3)*(edges(2,1)*edges(3,2) - edges(3,1)*edges(2,2)))/6.0_dp
        centroid = sum(v, dim=2)/4.0_dp
        if (coord_system.eq.1) then
            jacobian = centroid(1)
        else
            jacobian = abs(tetra_physics(i)%sqg1 + sum(tetra_physics(i)%gsqg*(centroid - v(:,1))))
        endif
        c%tetr_volume(i) = coordinate_volume*jacobian
    enddo

    print*, 'total physical volume of the tetrahedra = ', sum(c%tetr_volume)

end subroutine calc_physical_tetrahedron_volumes

subroutine build_reservoir_cells(min_volume_fraction)

    !Groups tetrahedra into reservoir cells for the energy and momentum conserving collisions. Every
    !tetrahedron is its own cell, except close to the magnetic axis: starting from the innermost
    !flux-surface layer and stopping at the first layer without tetrahedra smaller than
    !min_volume_fraction times the median tetrahedron volume, the tetrahedra of each layer and
    !toroidal slice are merged poloidally (in order of their poloidal angle) into cells of at least
    !that volume. Innermost layers that are too thin for this even as complete rings (their thinness
    !comes from the grid geometry near the axis) are first combined radially, from the axis outwards,
    !until the ring of every toroidal slice reaches the volume. A layer is identified by the lowest vertex flux of a tetrahedron, a toroidal slice by
    !its lowest vertex angle phi. min_volume_fraction <= 0 keeps one cell per tetrahedron.
    !All tetrahedra of a cell share the cell's flow and temperature (initialised to the
    !particle-weighted mean of its tetrahedra, exact while the flow is zero).

    use gorilla_applets_types_mod, only: c
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_rphiz, verts_sthetaphi
    use tetra_physics_mod, only: tetra_physics, coord_system
    use constants, only: pi

    real(dp), intent(in) :: min_volume_fraction

    real(dp), parameter :: tol = 1.0e-9_dp
    real(dp), dimension(:), allocatable :: s_low, phi_low, theta, centre_a, centre_b, n_group, sort_key, particles
    real(dp), dimension(:), allocatable :: a_c, b_c, ring_acc
    real(dp), dimension(:,:), allocatable :: ring_volume
    integer, dimension(:), allocatable :: radial_cell
    integer, dimension(:), allocatable :: layer_id, slice_id, group_id, perm, n_members
    logical, dimension(:), allocatable :: layer_has_small
    real(dp), dimension(3,4) :: v
    real(dp), dimension(4) :: s_v
    real(dp) :: threshold, cell_volume, u_mean, temp_mean
    integer :: i, k, j, idx, n_layers, n_slices, n_groups, n_cells, ind_phi, ind_a, ind_b, merge_limit, &
               group_first_cell, cell_first_idx, previous_group, n_merged_tetr, n_layers_grid
    logical :: radial_cell_open

    if (allocated(c%cell_of_tetr)) deallocate(c%cell_of_tetr)
    if (allocated(c%cell_start)) deallocate(c%cell_start)
    if (allocated(c%cell_members)) deallocate(c%cell_members)
    allocate(c%cell_of_tetr(ntetr))

    if (min_volume_fraction.le.0.0_dp) then
        c%cell_of_tetr = [(i, i = 1, ntetr)]
        n_cells = ntetr
    else
        ind_phi = 2
        ind_a = 1 !poloidal plane coordinates whose centroid angle orders the tetrahedra of a layer
        ind_b = 3
        if (coord_system.eq.2) ind_phi = 3
        allocate(s_low(ntetr), phi_low(ntetr), a_c(ntetr), b_c(ntetr), theta(ntetr))
        do i = 1, ntetr
            do k = 1, 4
                if (coord_system.eq.1) then
                    v(:,k) = verts_rphiz(:,tetra_grid(i)%ind_knot(k))
                else
                    v(:,k) = verts_sthetaphi(:,tetra_grid(i)%ind_knot(k))
                endif
            enddo
            v(ind_phi,2:4) = v(ind_phi,1) + modulo(v(ind_phi,2:4) - v(ind_phi,1) + pi, 2.0_dp*pi) - pi
            if (coord_system.eq.1) then
                do k = 1, 4
                    s_v(k) = normalised_poloidal_flux(tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi*(v(:,k) - v(:,1))))
                enddo
                a_c(i) = sum(v(ind_a,:))/4.0_dp
                b_c(i) = sum(v(ind_b,:))/4.0_dp
            else
                s_v = v(1,:)
                !flux coordinates: order by the poloidal angle theta itself
                v(2,2:4) = v(2,1) + modulo(v(2,2:4) - v(2,1) + pi, 2.0_dp*pi) - pi
                a_c(i) = cos(sum(v(2,:))/4.0_dp)
                b_c(i) = sin(sum(v(2,:))/4.0_dp)
            endif
            s_low(i) = minval(s_v)
            phi_low(i) = modulo(minval(v(ind_phi,:)), 2.0_dp*pi)
            if (phi_low(i).gt.2.0_dp*pi - tol) phi_low(i) = 0.0_dp
        enddo

        allocate(layer_id(ntetr), slice_id(ntetr))
        call cluster_values(s_low, tol, layer_id, n_layers)
        call cluster_values(phi_low, tol, slice_id, n_slices)
        allocate(group_id(ntetr))
        threshold = min_volume_fraction*median(c%tetr_volume)

        !radial merging of layers whose rings are too small, from the axis outwards
        n_layers_grid = n_layers
        allocate(ring_volume(n_layers, n_slices), ring_acc(n_slices), radial_cell(n_layers))
        ring_volume = 0.0_dp
        do i = 1, ntetr
            ring_volume(layer_id(i), slice_id(i)) = ring_volume(layer_id(i), slice_id(i)) + c%tetr_volume(i)
        enddo
        n_layers = 0
        radial_cell_open = .false.
        do k = 1, n_layers_grid
            if (.not.radial_cell_open) then
                n_layers = n_layers + 1
                ring_acc = 0.0_dp
                radial_cell_open = .true.
            endif
            radial_cell(k) = n_layers
            ring_acc = ring_acc + ring_volume(k,:)
            if (minval(ring_acc).ge.threshold) radial_cell_open = .false.
        enddo
        layer_id = radial_cell(layer_id)

        group_id = (layer_id - 1)*n_slices + slice_id
        n_groups = n_layers*n_slices

        !poloidal angle of each centroid around the centre of its layer and slice
        allocate(centre_a(n_groups), centre_b(n_groups), n_group(n_groups))
        centre_a = 0.0_dp
        centre_b = 0.0_dp
        n_group = 0.0_dp
        do i = 1, ntetr
            centre_a(group_id(i)) = centre_a(group_id(i)) + a_c(i)
            centre_b(group_id(i)) = centre_b(group_id(i)) + b_c(i)
            n_group(group_id(i)) = n_group(group_id(i)) + 1.0_dp
        enddo
        centre_a = centre_a/max(n_group, 1.0_dp)
        centre_b = centre_b/max(n_group, 1.0_dp)
        if (coord_system.eq.2) then
            centre_a = 0.0_dp
            centre_b = 0.0_dp
        endif
        do i = 1, ntetr
            theta(i) = atan2(b_c(i) - centre_b(group_id(i)), a_c(i) - centre_a(group_id(i)))
        enddo

        !layers 1 .. merge_limit (from the axis outwards) contain small tetrahedra and are merged
        allocate(layer_has_small(n_layers))
        layer_has_small = .false.
        do i = 1, ntetr
            if (c%tetr_volume(i).lt.threshold) layer_has_small(layer_id(i)) = .true.
        enddo
        merge_limit = 0
        do while (merge_limit.lt.n_layers)
            if (.not.layer_has_small(merge_limit+1)) exit
            merge_limit = merge_limit + 1
        enddo

        allocate(sort_key(ntetr), perm(ntetr))
        sort_key = real(group_id, dp) + 0.999_dp*(theta + pi)/(2.0_dp*pi)
        call sort_index(sort_key, perm)

        n_cells = 0
        previous_group = -1
        group_first_cell = 0
        cell_first_idx = 1
        cell_volume = 0.0_dp
        do idx = 1, ntetr
            i = perm(idx)
            if (layer_id(i).gt.merge_limit) then
                n_cells = n_cells + 1
                c%cell_of_tetr(i) = n_cells
                previous_group = group_id(i)
                cycle
            endif
            if ((group_id(i).ne.previous_group).or.(cell_volume.ge.threshold)) then
                n_cells = n_cells + 1
                if (group_id(i).ne.previous_group) group_first_cell = n_cells
                cell_first_idx = idx
                cell_volume = 0.0_dp
            endif
            c%cell_of_tetr(i) = n_cells
            cell_volume = cell_volume + c%tetr_volume(i)
            previous_group = group_id(i)
            !at the end of a merged group, a last cell below the threshold joins the previous one
            if (idx.lt.ntetr) then
                if (group_id(perm(idx+1)).eq.group_id(i)) cycle
            endif
            if ((cell_volume.lt.threshold).and.(n_cells.gt.group_first_cell)) then
                c%cell_of_tetr(perm(cell_first_idx:idx)) = n_cells - 1
                n_cells = n_cells - 1
            endif
            cell_volume = threshold !the next tetrahedron starts a new cell
        enddo
    endif

    !cell members in compressed form
    allocate(n_members(n_cells), c%cell_start(n_cells+1), c%cell_members(ntetr))
    n_members = 0
    do i = 1, ntetr
        n_members(c%cell_of_tetr(i)) = n_members(c%cell_of_tetr(i)) + 1
    enddo
    c%cell_start(1) = 1
    do k = 1, n_cells
        c%cell_start(k+1) = c%cell_start(k) + n_members(k)
    enddo
    n_members = 0
    do i = 1, ntetr
        k = c%cell_of_tetr(i)
        c%cell_members(c%cell_start(k) + n_members(k)) = i
        n_members(k) = n_members(k) + 1
    enddo

    !common state of the merged cells
    n_merged_tetr = 0
    do k = 1, n_cells
        if (c%cell_start(k+1) - c%cell_start(k).le.1) cycle
        n_merged_tetr = n_merged_tetr + c%cell_start(k+1) - c%cell_start(k)
        do j = 1, c%n
            associate(members => c%cell_members(c%cell_start(k):c%cell_start(k+1)-1))
                particles = c%dens_mat(j,members)*c%tetr_volume(members)
                u_mean = sum(particles*c%vpar_mat(j,members))/sum(particles)
                temp_mean = sum(particles*c%temp_mat(j,members))/sum(particles)
                c%vpar_mat(j,members) = u_mean
                c%temp_mat(j,members) = temp_mean
            end associate
        enddo
    enddo

    if (min_volume_fraction.gt.0.0_dp) then
        print*, 'reservoir cells: ', n_cells, ' cells for ', ntetr, ' tetrahedra; ', n_merged_tetr, &
                ' tetrahedra in ', merge_limit, ' layers merged poloidally, threshold volume ', threshold
        print*, 'innermost ', n_layers_grid - n_layers, ' grid layers merged radially (', &
                n_layers_grid, ' grid layers ->', n_layers, ')'
        print*, 'smallest reservoir cell volume = ', minval([(sum(c%tetr_volume(c%cell_members(c%cell_start(k): &
                c%cell_start(k+1)-1))), k = 1, n_cells)])
    endif

contains

    function median(values) result(m)
        real(dp), dimension(:), intent(in) :: values
        real(dp) :: m
        integer, dimension(size(values)) :: p
        call sort_index(values, p)
        m = values(p((size(values)+1)/2))
    end function median

end subroutine build_reservoir_cells

subroutine cluster_values(values, tol, ids, n_ids)

    !ids(i) = rank of the cluster of values(i) among clusters of values closer than tol, ascending

    real(dp), dimension(:), intent(in) :: values
    real(dp), intent(in) :: tol
    integer, dimension(:), intent(out) :: ids
    integer, intent(out) :: n_ids
    integer, dimension(size(values)) :: p
    integer :: k

    call sort_index(values, p)
    n_ids = 1
    ids(p(1)) = 1
    do k = 2, size(values)
        if (values(p(k)) - values(p(k-1)).gt.tol) n_ids = n_ids + 1
        ids(p(k)) = n_ids
    enddo

end subroutine cluster_values

subroutine sort_index(a, p)

    !Heapsort: p is the permutation that sorts a ascending

    real(dp), dimension(:), intent(in) :: a
    integer, dimension(:), intent(out) :: p
    integer :: n, i, last, tmp

    n = size(a)
    p = [(i, i = 1, n)]
    do i = n/2, 1, -1
        call sift_down(i, n)
    enddo
    do last = n, 2, -1
        tmp = p(1)
        p(1) = p(last)
        p(last) = tmp
        call sift_down(1, last-1)
    enddo

contains

    subroutine sift_down(start, end)
        integer, intent(in) :: start, end
        integer :: root, child, t
        root = start
        do while (2*root.le.end)
            child = 2*root
            if (child.lt.end) then
                if (a(p(child)).lt.a(p(child+1))) child = child + 1
            endif
            if (a(p(root)).ge.a(p(child))) return
            t = p(root)
            p(root) = p(child)
            p(child) = t
            root = child
        enddo
    end subroutine sift_down

end subroutine sort_index

subroutine deallocate_collision_arrays
!
! Deallocates collision coefficient arrays to allow reallocation when grid parameters change.
!
    use gorilla_applets_types_mod, only: c

    if (allocated(c%dens_mat))   deallocate(c%dens_mat)
    if (allocated(c%temp_mat))   deallocate(c%temp_mat)
    if (allocated(c%vpar_mat))   deallocate(c%vpar_mat)
    if (allocated(c%efcolf_mat)) deallocate(c%efcolf_mat)
    if (allocated(c%velrat_mat)) deallocate(c%velrat_mat)
    if (allocated(c%enrat_mat))  deallocate(c%enrat_mat)
    if (allocated(c%mass))       deallocate(c%mass)
    if (allocated(c%charge_num)) deallocate(c%charge_num)
    if (allocated(c%dens))       deallocate(c%dens)
    if (allocated(c%temp))       deallocate(c%temp)
    if (allocated(c%randcol))    deallocate(c%randcol)
    if (allocated(c%tetr_volume)) deallocate(c%tetr_volume)
    if (allocated(c%cell_of_tetr)) deallocate(c%cell_of_tetr)
    if (allocated(c%cell_start)) deallocate(c%cell_start)
    if (allocated(c%cell_members)) deallocate(c%cell_members)

end subroutine deallocate_collision_arrays

subroutine set_c

    use gorilla_applets_types_mod, only: c, in
    use tetra_grid_mod, only: ntetr
    use constants, only: amp, ame

    c%n = in%n_background_species

    !Free any pre-existing collision arrays so a fresh build can proceed without external bookkeeping.
    call deallocate_collision_arrays()

    allocate(c%dens_mat(c%n, ntetr))
    allocate(c%temp_mat(c%n, ntetr))
    allocate(c%vpar_mat(c%n, ntetr))
    allocate(c%efcolf_mat(c%n, ntetr))
    allocate(c%velrat_mat(c%n, ntetr))
    allocate(c%enrat_mat(c%n, ntetr))
    allocate(c%mass(c%n))
    allocate(c%charge_num(c%n))
    allocate(c%dens(c%n))
    allocate(c%temp(c%n))

    c%mass(1:c%n-1) = 2.0_dp * amp
    c%mass(c%n) = ame
    c%charge_num(1:c%n-1) = 1.0_dp
    c%charge_num(c%n) = -1.0_dp
    c%vpar_mat = 0.0_dp

end subroutine set_c

subroutine get_T_and_n_from_files

    use gorilla_applets_types_mod, only: c
    use tetra_grid_settings_mod, only: grid_size
    use tetra_grid_mod, only: ntetr

    integer :: Te_unit, Ti_unit, ne_unit, i
    
    open(newunit = Te_unit, file = 'background/Te_d.dat')
    read(Te_unit,'(e16.9)') (c%temp_mat(c%n, i), i=1, ntetr/grid_size(2), 3)
    close(Te_unit)

    open(newunit = Ti_unit, file = 'background/Ti_d.dat')
    read(Ti_unit,'(e16.9)') (c%temp_mat(1, i), i=1, ntetr/grid_size(2), 3)
    close(Ti_unit)

    open(newunit = ne_unit, file = 'background/ne_d.dat')
    read(ne_unit,'(e16.9)') (c%dens_mat(1, i), i=1, ntetr/grid_size(2), 3)
    c%dens_mat(c%n, :) = c%dens_mat(1, :)
    close(ne_unit)

end subroutine get_T_and_n_from_files

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

subroutine normalise_prism_moments_and_prism_moments_squared(species_in, boole_skip_time_normalisation_in)

    use gorilla_applets_types_mod, only: moment_specs, output, in, start

    integer, intent(in), optional :: species_in
    logical, intent(in), optional :: boole_skip_time_normalisation_in
    integer :: species
    integer :: n
    real(dp) :: time
    logical :: boole_skip_time_normalisation

    time = in%time_step
    species = 1
    if (present(species_in)) species = species_in
    time = start%t(species)

    boole_skip_time_normalisation = .false.
    if (present(boole_skip_time_normalisation_in)) boole_skip_time_normalisation = boole_skip_time_normalisation_in

    if (boole_skip_time_normalisation) time = 1.0_dp

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
    complex(dp)                                 :: i
    complex(dp), dimension(:,:,:), allocatable  :: prism_moments_ordered_for_ft
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

    use gorilla_applets_types_mod, only: in, weights

    integer  :: i
    real(dp) :: maximum_weight, minimum_weight, average_weight

    maximum_weight = abs(weights%w(1,1))
    minimum_weight = abs(weights%w(1,1))
    average_weight = 0

    do i = 1,int(in%n_particles)
        average_weight = average_weight + abs(weights%w(i,1))
        maximum_weight = max(maximum_weight, abs(weights%w(i,1)))
        minimum_weight = min(minimum_weight, abs(weights%w(i,1)))
    enddo

    average_weight = average_weight/in%n_particles

    print*, "maximum particle weight = ", maximum_weight
    print*, "minimum particle weight = ", minimum_weight
    print*, "average particle weight = ", average_weight

end subroutine analyse_particle_weight_distribution

end module utils_data_pre_and_post_processing_mod