module boltzmann_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine read_boltzmann_inp_into_type

    use boltzmann_types_mod, only: in

    real(dp) :: time_step,energy_eV,n_particles, density
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_boltzmann_energies, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
               boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
               boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions
    integer :: i_integrator_type, seed_option, n_background_density_updates

    integer :: b_inp_unit

    !Namelist for boltzmann input
    NAMELIST /boltzmann_nml/ time_step,energy_eV,n_particles,boole_squared_moments,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, &
    & boole_antithetic_variate,boole_linear_temperature_simulation,i_integrator_type,seed_option, boole_write_vertex_indices, &
    & boole_write_vertex_coordinates, boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
    & boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_background_density_updates

    open(newunit = b_inp_unit, file='boltzmann.inp', status='unknown')
    read(b_inp_unit,nml=boltzmann_nml)
    close(b_inp_unit)

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
    in%n_background_density_updates = n_background_density_updates

    print *,'GORILLA: Loaded input data from boltzmann.inp'

end subroutine read_boltzmann_inp_into_type

subroutine calc_boltzmann

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg,pi,echarge
    use tetra_physics_mod, only: particle_mass,particle_charge
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_grid_mod, only: ntetr
    use gorilla_settings_mod, only: ispecies
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use boltzmann_types_mod, only: moment_specs, counter, c, in
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, &
    get_ipert, set_moment_specifications, initialise_output, calc_starting_conditions, initialize_exit_data, calc_poloidal_flux, &
    calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared, fourier_transform_moments, &
    find_minimal_angle_between_curlA_and_tetrahedron_faces, analyse_particle_weight_distribution, &
    perform_background_density_update, set_weights
    use boltzmann_types_mod, only: output

    real(dp) :: v0
    real(dp), dimension(:,:), allocatable :: verts
    integer :: i

    call set_seed_for_random_numbers
    call read_boltzmann_inp_into_type
    call get_ipert()
    call initialize_gorilla(i_option,ipert)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete this again afterwards !!!!!!!!!!!!!!!!!!!!!!!
    if (ispecies.eq.4) particle_charge = 15*echarge
    print*, 'particle charge number = ', particle_charge/echarge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    v0=sqrt(2.0_dp*in%energy_eV*ev2erg/particle_mass)

    call set_moment_specifications
    call initialise_output
    call calc_square_root_g
    call calc_volume_integrals(in%boole_boltzmann_energies,in%boole_refined_sqrt_g, in%density, in%energy_eV)
    call initialize_exit_data
    call calc_starting_conditions(verts)
    call calc_poloidal_flux(verts)
    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra(v0)
    call give_file_names
    call unlink_files

    if (in%i_integrator_type.eq.2) print*, 'Error: i_integratpr_type set to 2, this module only works with &
                                    & i_integrator_type set to 1'

    if (in%n_background_density_updates.eq.0) then
        call parallelised_particle_pushing(v0)
    else
        do i = 1, in%n_background_density_updates
            !set everything 0, also counters, write subroutine for that
            !potentially set weights to updated densities
            call set_weights !weights need to be set again because in orbit_timestep_gorilla_boltzmannn they are multiplied with 
            !some factor, so we need to get rid of it again
            output%tetr_moments = 0.0_dp
            output%prism_moments = 0.0_dp
            output%prism_moments_squared = 0.0_dp
            call parallelised_particle_pushing(v0)
            call perform_background_density_update(i)
        enddo
    endif

    call normalise_prism_moments_and_prism_moments_squared
    if (moment_specs%n_moments.gt.0) call fourier_transform_moments
    call write_data_to_files

    if (in%boole_precalc_collisions) print*, "maxcol = ", c%maxcol
    print*, 'Number of lost particles',counter%lost_particles
    print*, 'average number of pushings = ', counter%tetr_pushings/in%n_particles
    print*, 'average number of toroidal revolutions = ', counter%phi_0_mappings/in%n_particles
    print*, 'average number of integration steps = ', counter%integration_steps/in%n_particles
    PRINT*, 'particle mass = ', particle_mass
    PRINT*, 'absolute value of velocity = ', v0
    PRINT*, 'particle charge = ', particle_charge
    PRINT*, 'temperature = ', ev2erg*in%energy_eV
    print*, 'energy in eV = ', in%energy_eV
    print*, 'tracing time in seconds = ', in%time_step
    if((grid_kind.eq.2).or.(grid_kind.eq.3)) then
         print*, 'number of times particles were pushed across the inside hole = ', counter%lost_inside
    endif

end subroutine calc_boltzmann

subroutine parallelised_particle_pushing(v0)

    use boltzmann_types_mod, only: counter, c, in, time_t, moment_specs, counter_t, particle_status_t
    use tetra_grid_mod, only: ntetr
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use utils_parallelised_particle_pushing_mod, only: print_progress, handle_lost_particles, add_local_tetr_moments_to_output, &
    add_local_counter_to_counter, initialise_loop_variables, carry_out_collisions, update_exit_data, &
    initialise_seed_for_random_numbers_for_each_thread

    real(dp), intent(in)                     :: v0
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
    !$OMP& SHARED(counter, kpart,v0, in, c, iantithetic) &
    !$OMP& PRIVATE(p,l,n,i,x,vpar,vperp,t,ind_tetr,iface,local_tetr_moments,local_counter,particle_status) &
    !$OMP& FIRSTPRIVATE(thread_flag)
    print*, 'get number of threads', omp_get_num_threads()
    !$OMP DO SCHEDULE(static)

    !Loop over particles
    do p = 1,in%num_particles/iantithetic

        if ((.not.in%boole_precalc_collisions).and.thread_flag) then
            call initialise_seed_for_random_numbers_for_each_thread(omp_get_thread_num())
            thread_flag = .false.
        endif

        do l = 1,iantithetic
            n = (p-1)*iantithetic+l
            !$omp critical
            kpart = kpart+1 !in general not equal to n becuase of parallelisation
            call print_progress(in%num_particles,kpart,n)
            !$omp end critical

            call initialise_loop_variables(l, n, v0, local_counter,particle_status,t,local_tetr_moments,x,vpar,vperp)

            i = 0
            do while (t%confined.lt.in%time_step)
                i = i+1

                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, v0, t, x, vpar,vperp,ind_tetr, iface)
                    t%step = t%step/v0 !in carry_out_collisions, t%step is initiated as a length, so you need to divide by v0
                endif

                call orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t%step,particle_status,ind_tetr,iface,n,&
                            & local_tetr_moments, local_counter,t%remain)

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
            call update_exit_data(particle_status%lost,t%confined,x,vpar,vperp,i,n)
        enddo
        !$omp critical
        call add_local_tetr_moments_to_output(local_tetr_moments)
        !$omp end critical
    enddo !n
    !$OMP END DO
    !$OMP END PARALLEL

end subroutine

subroutine orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,particle_status,ind_tetr,iface, n,local_tetr_moments, &
                                            local_counter, t_remain_out)

    use pusher_tetra_rk_mod, only: pusher_tetra_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: vperp_func
    use find_tetra_mod, only: find_tetra
    use boltzmann_types_mod, only: counter_t, particle_status_t, g
    use tetra_grid_settings_mod, only: grid_kind
    use utils_orbit_timestep_mod, only: identify_particles_entering_annulus, update_local_tetr_moments, &
                                                       initialize_constants_of_motion, calc_particle_weights_and_jperp

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
        call calc_particle_weights_and_jperp(n,z_save,vpar,vperp,ind_tetr)
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
        if(boole_t_finished) then !Orbit stops within cell, because "flight"-time t_step has finished
            if( present(t_remain_out)) t_remain_out = t_remain
            exit
        endif

    enddo !Loop for tetrahedron pushings

    vperp = vperp_func(z_save,perpinv,ind_tetr_save) !Compute vperp from position

end subroutine orbit_timestep_gorilla_boltzmann

end module boltzmann_mod