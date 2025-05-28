module self_consistent_electric_field_mod

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
               boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions
    integer :: i_integrator_type, seed_option, n_electric_potential_updates, update_dimension, n_species

    integer :: s_inp_unit

    !Namelist for self consistent electric field input
    NAMELIST /self_consistent_ef_nml/ time_step,energy_eV,n_particles,boole_squared_moments,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, &
    & boole_antithetic_variate,boole_linear_temperature_simulation,i_integrator_type,seed_option, boole_write_vertex_indices, &
    & boole_write_vertex_coordinates, boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
    & boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_electric_potential_updates, update_dimension, &
    & n_species

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

    print *,'GORILLA_APPLETS: Loaded input data from self_consistent_ef.inp'

end subroutine read_self_consistent_electric_field_inp_into_type

subroutine calc_self_consistent_electric_field

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg,pi,echarge
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_grid_mod, only: ntetr, verts_sthetaphi
    use gorilla_settings_mod, only: ispecies
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_volume_integrals_in_flux_coordinates
    use gorilla_applets_types_mod, only: moment_specs, counter, c, in, start
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, &
    get_ipert, set_moment_specifications, initialise_output, initialize_exit_data, calc_poloidal_flux, &
    calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared, fourier_transform_moments, &
    find_minimal_angle_between_curlA_and_tetrahedron_faces, analyse_particle_weight_distribution, &
    perform_electric_potential_update, set_weights, prepare_next_round_of_parallelised_particle_pushing, &
    associate_flux_labels_with_tetrahedra_and_vertices, allocate_electric_potential_type
    use gorilla_applets_types_mod, only: output, ep

    integer :: i, species

    call set_seed_for_random_numbers
    call read_self_consistent_electric_field_inp_into_type
    call get_ipert()
    call initialize_gorilla(i_option,ipert)
    call associate_flux_labels_with_tetrahedra_and_vertices
    call set_moment_specifications
    call initialise_output
    call calc_volume_integrals_in_flux_coordinates
    call initialize_exit_data
    call calc_starting_conditions
    call calc_poloidal_flux(verts_sthetaphi)
    call allocate_electric_potential_type
    call give_file_names
    call unlink_files
    call print_errors_for_bad_inputs

    do i = 1, max(in%n_electric_potential_updates,1)
        ep%rho_prism = 0
        do species = 1,2 !trace electrons and ions
            call prepare_next_round_of_parallelised_particle_pushing(species)
            if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra(species)
            call parallelised_particle_pushing(species,i)
            call normalise_prism_moments_and_prism_moments_squared(species)
            ep%rho_prism = ep%rho_prism + real(output%prism_moments(1,:,species))*start%particle_charge(species)
        enddo
        call perform_electric_potential_update(i)
    enddo

    if (moment_specs%n_moments.gt.0) call fourier_transform_moments
    call write_data_to_files

    if (in%boole_precalc_collisions) print*, "maxcol = ", c%maxcol
    print*, 'Number of lost ions',counter%lost_particles
    print*, 'average number of pushings = ', counter%tetr_pushings/in%n_particles
    print*, 'average number of toroidal revolutions = ', counter%phi_0_mappings/in%n_particles
    print*, 'average number of integration steps = ', counter%integration_steps/in%n_particles
    PRINT*, 'ion mass = ', start%particle_mass(1)
    PRINT*, 'absolute value of ion velocity = ', start%v0(1)
    PRINT*, 'ion charge = ', start%particle_charge(1)
    PRINT*, 'temperature = ', ev2erg*in%energy_eV
    print*, 'energy in eV = ', in%energy_eV
    print*, 'tracing time in seconds = ', in%time_step
    if((grid_kind.eq.2).or.(grid_kind.eq.3)) then
         print*, 'number of times that particles were pushed across the inside hole = ', counter%lost_inside
    endif
    print*, 'Average abs Delta Phi at all the electric potential updates = ', ep%average_abs_phi_elec_from_rho

end subroutine calc_self_consistent_electric_field

subroutine parallelised_particle_pushing(species,j)

    use gorilla_applets_types_mod, only: counter, c, in, time_t, moment_specs, counter_t, particle_status_t, start, exit_data
    use tetra_grid_mod, only: ntetr
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use utils_parallelised_particle_pushing_mod, only: print_progress, handle_lost_particles, add_local_tetr_moments_to_output, &
    add_local_counter_to_counter, initialise_loop_variables, carry_out_collisions, update_exit_data, update_start_type, &
    initialise_seed_for_random_numbers_for_each_thread
    

    integer, intent(in)                      :: species, j
    integer                                  :: kpart, iantithetic, ind_tetr, iface
    integer                                  :: p, l, n, i
    real(dp), dimension(3)                   :: x
    real(dp)                                 :: vpar,vperp
    type(time_t)                             :: t
    type(counter_t)                          :: local_counter
    type(particle_status_t)                  :: particle_status
    complex(dp), dimension(:,:), allocatable :: local_tetr_moments
    logical                                  :: thread_flag = .true.

    allocate(local_tetr_moments(moment_specs%n_moments,ntetr))
    kpart = 0
    iantithetic = 1
    if (in%boole_antithetic_variate) iantithetic = 2

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(counter, kpart,species, in, c, iantithetic, start, j) &
    !$OMP& PRIVATE(p,l,n,i,x,vpar,vperp,t,ind_tetr,iface,local_tetr_moments,local_counter,particle_status) &
    !$OMP& FIRSTPRIVATE(thread_flag)
    print*, 'get number of threads', omp_get_num_threads()
    !$OMP DO SCHEDULE(static)

    !Loop over particles
    do p = 1,in%num_particles/iantithetic

        if ((.not.in%boole_precalc_collisions).and.thread_flag) then
            call initialise_seed_for_random_numbers_for_each_thread(omp_get_thread_num(), j)
            thread_flag = .false.
        endif

        do l = 1,iantithetic
            n = (p-1)*iantithetic+l
            !$omp critical
            kpart = kpart+1 !in general not equal to n becuase of parallelisation
            call print_progress(in%num_particles,kpart,n)
            !$omp end critical

            call initialise_loop_variables(l, n, local_counter,particle_status,t,local_tetr_moments,x,vpar,vperp,species)

            i = 0
            do while (t%confined.lt.start%t(species))
                i = i+1
                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, t, x, vpar,vperp,ind_tetr, iface, species)
                    t%step = t%step/start%v0(species) !in carry_out_collisions, t%step is initiated as a length, so you need to divide by v0
                endif
                call orbit_timestep_gorilla_self_consistent_ef(x,vpar,vperp,t%step,particle_status,ind_tetr,iface,n,&
                            & local_tetr_moments, local_counter,t%remain, species, j)

                t%confined = t%confined + t%step - t%remain

                if (ind_tetr.eq.-1) then
                    call handle_lost_particles(local_counter, particle_status%lost)
                    exit
                endif
            enddo


            !$omp critical
            counter%integration_steps = counter%integration_steps + i
            c%maxcol = max(dble(i)/dble(c%randcoli),c%maxcol)
            call add_local_counter_to_counter(local_counter)
            !$omp end critical
            call update_exit_data(particle_status%lost,t%confined,x,vpar,vperp,i,n,species_in=species)
            call update_start_type(x,vpar,vperp,n,species,ind_tetr)
        enddo
        !$omp critical
        call add_local_tetr_moments_to_output(local_tetr_moments,species)
        !$omp end critical
    enddo !n
    !$OMP END DO
    !$OMP END PARALLEL

end subroutine parallelised_particle_pushing

subroutine orbit_timestep_gorilla_self_consistent_ef(x,vpar,vperp,t_step,particle_status,ind_tetr,iface, n,local_tetr_moments, &
                                            local_counter, t_remain_out,species, j)

    use pusher_tetra_rk_mod, only: pusher_tetra_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: vperp_func
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_types_mod, only: counter_t, particle_status_t, g, start
    use tetra_grid_settings_mod, only: grid_kind
    use utils_orbit_timestep_mod, only: identify_particles_entering_annulus, update_local_tetr_moments, &
                                        initialize_constants_of_motion, compute_radial_fluxes

    integer, intent(in)                          :: species, j
    type(counter_t), intent(inout)               :: local_counter
    type(particle_status_t), intent(inout)       :: particle_status
    real(dp), dimension(3), intent(inout)        :: x
    complex(dp), dimension(:,:), intent (inout)  :: local_tetr_moments
    real(dp), intent(inout)                      :: vpar,vperp
    real(dp), intent(in)                         :: t_step
    integer, intent(inout)                       :: ind_tetr,iface
    real(dp), intent(out), optional              :: t_remain_out
    real(dp), dimension(3)                       :: z_save, x_save
    real(dp)                                     :: t_remain,t_pass,perpinv
    logical                                      :: boole_t_finished, boole_lost_inside
    integer                                      :: ind_tetr_save,iper_phi,n
    type(optional_quantities_type)               :: optional_quantities
    real(dp), dimension(3)                       :: x_new
    
    if(.not.particle_status%initialized) then !If orbit_timestep is called for the first time without grid position
        call check_coordinate_domain(x) !Check coordinate domain (optionally perform modulo operation)
        call find_tetra(x,vpar,vperp,ind_tetr,iface) !Find tetrahedron index and face index for position x
        if(ind_tetr.eq.-1) then !If particle doesn't lie inside any tetrahedron
            t_remain_out = t_step
            return
        endif
        z_save = x-tetra_physics(ind_tetr)%x1
        if (j.eq.1) call calc_particle_weights_and_jperp(n,z_save,vpar,vperp,ind_tetr,species)
        particle_status%initialized = .true.
    endif
          
    if(t_step.eq.0.0_dp) return !Exit the subroutine after initialization, if time step equals zero
    if(particle_status%initialized) z_save = x-tetra_physics(ind_tetr)%x1
    call initialize_constants_of_motion(vperp,z_save,ind_tetr,perpinv)
    t_remain = t_step
    boole_t_finished = .false.
    local_counter%tetr_pushings = local_counter%tetr_pushings -1 !set tetr_pushings to -1 because when entering the loop it will go back to one without pushing

    do !Loop for tetrahedron pushings until t_step is reached
        local_counter%tetr_pushings = local_counter%tetr_pushings +1

        if(ind_tetr.eq.-1) then
            if(present(t_remain_out)) t_remain_out = t_remain
            if((grid_kind.eq.2).or.(grid_kind.eq.3)) then
                call identify_particles_entering_annulus(x,local_counter,boole_lost_inside)
                if (boole_lost_inside) then
                    x_new = 3*(/g%raxis,x(2),g%zaxis/) - 2*x
                    call find_tetra(x_new,vpar,vperp,ind_tetr,iface)
                    x = x_new
                    print*, "particle pushing across the hole surrounding the magnetic axis was successful"
                    if (ind_tetr.eq.-1) then
                        print*, "ATTENTION: particle pushing across the hole surrounding the magnetic axis was unsuccessful"
                        exit
                    endif
                else
                    exit
                endif
            else
                exit
            endif
        endif

        ind_tetr_save = ind_tetr
        x_save = x

        select case(ipusher) !Calculate trajectory
            case(1)
                call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper_phi)
            case(2) 
                call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_remain,&
                                                    & t_pass,boole_t_finished,iper_phi,optional_quantities)
        end select

        t_remain = t_remain - t_pass

        call update_local_tetr_moments(local_tetr_moments,ind_tetr_save,n,optional_quantities)
        if((grid_kind.eq.2).or.(grid_kind.eq.3)) call compute_radial_fluxes(ind_tetr_save,ind_tetr,x)

        if(boole_t_finished) then !Orbit stops within cell, because "flight"-time t_step has finished
            if( present(t_remain_out)) t_remain_out = t_remain
            exit
        endif

    enddo !Loop for tetrahedron pushings

    vperp = vperp_func(z_save,perpinv,ind_tetr_save) !Compute vperp from position

end subroutine orbit_timestep_gorilla_self_consistent_ef

subroutine print_errors_for_bad_inputs

    use gorilla_applets_types_mod, only: in

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

end subroutine print_errors_for_bad_inputs

subroutine calc_particle_weights_and_jperp(n,z_save,vpar,vperp,ind_tetr, species_in)

    use gorilla_applets_types_mod, only: in, flux, start
    use tetra_physics_mod, only: tetra_physics
    use constants, only: ev2erg, pi
    use volume_integrals_and_sqrt_g_mod, only: sqrt_g
    use supporting_functions_mod, only: bmod_func

    real(dp), intent(in) :: vpar, vperp
    real(dp), dimension(3), intent(in) :: z_save
    integer, intent(in) :: n,ind_tetr
    integer, intent(in), optional :: species_in
    integer :: species = 1
    real(dp) :: local_poloidal_flux, phi_elec_func, temperature, epsilon_max

    if (present(species_in)) species = species_in

    start%weight(n,species) = start%weight(n,species)*abs((tetra_physics(ind_tetr)%sqg1 + sum(tetra_physics(ind_tetr)%gsqg*z_save)))

    if (in%boole_linear_density_simulation.or.in%boole_linear_temperature_simulation) then
        local_poloidal_flux = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*z_save)
    endif
    if (in%boole_linear_density_simulation) then
        start%weight(n,species) = start%weight(n,species)*(flux%poloidal_max*1.1_dp-local_poloidal_flux)/(flux%poloidal_max*1.1_dp)
    endif

    if (in%boole_boltzmann_energies) then
        phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z_save)
        epsilon_max = 5*in%energy_eV*ev2erg
        start%weight(n,species) = start%weight(n,species)*epsilon_max*2/sqrt(pi)*sqrt(start%energy(n,species)*ev2erg)
        if (.not. in%boole_linear_temperature_simulation) then
            start%weight(n,species) =start%weight(n,species) /(in%energy_eV*ev2erg)**1.5_dp* &
                        & exp(-(start%energy(n,species)*ev2erg+start%particle_charge(species)*phi_elec_func)/(in%energy_eV*ev2erg))
        else
            temperature = in%energy_eV*ev2erg*(flux%poloidal_max*1.1_dp-local_poloidal_flux)/(flux%poloidal_max*1.1_dp)
            start%weight(n,species) = start%weight(n,species)/temperature**1.5_dp* &
            & exp(-(start%energy(n,species)*ev2erg+start%particle_charge(species)*phi_elec_func)/temperature)
        endif
    endif

    start%jperp(n,species) = start%particle_mass(species)*vperp**2*start%cm_over_e(species)/(2*bmod_func(z_save,ind_tetr))*(-1)
    !-1 because of negative gyrophase

end subroutine calc_particle_weights_and_jperp

subroutine calc_starting_conditions

    use gorilla_applets_types_mod, only: in
    
    real(dp), dimension(:,:,:), allocatable                :: rand_matrix

    allocate(rand_matrix(5,in%num_particles,in%n_species))
    call RANDOM_NUMBER(rand_matrix)

    call allocate_start_type
    call set_starting_positions(rand_matrix)
    call set_rest_of_start_type(rand_matrix)

end subroutine calc_starting_conditions

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
    use tetra_grid_settings_mod, only: sfc_s_min, n_field_periods

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix

    start%x(1,:,:) = sfc_s_min + rand_matrix(1,:,:)*(1-sfc_s_min) !s
    start%x(2,:,:) = 2*pi*rand_matrix(2,:,:) !theta
    start%x(3,:,:) = 2*pi/n_field_periods*rand_matrix(3,:,:) !phi

end subroutine set_starting_positions

subroutine set_rest_of_start_type(rand_matrix)

    use gorilla_applets_types_mod, only: in, start
    use tetra_physics_mod, only: cm_over_e, particle_charge, particle_mass
    use gorilla_applets_settings_mod, only: i_option
    use constants, only: echarge,ame,clight, ev2erg, pi
    use tetra_grid_settings_mod, only: sfc_s_min, n_field_periods

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

    start%particle_charge = (/particle_charge, -echarge/)
    start%particle_mass = (/particle_mass, ame/)
    start%cm_over_e = (/cm_over_e, -clight*ame/echarge/)
    start%t = (/in%time_step, in%time_step/42.0_dp/)

    start%v0 = sqrt(2.0_dp*in%energy_eV*ev2erg/start%particle_mass)

    start%weight = in%density*(1-sfc_s_min)*4*pi**2/n_field_periods

end subroutine set_rest_of_start_type

end module self_consistent_electric_field_mod