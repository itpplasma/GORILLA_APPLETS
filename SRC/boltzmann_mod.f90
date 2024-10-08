module boltzmann_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine read_boltzmann_inp_into_type(b)

    use boltzmann_types_mod, only: boltzmann_input_t

    real(dp) :: time_step,energy_eV,n_particles, density
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_boltzmann_energies, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation
    integer :: i_integrator_type, seed_option

    integer :: b_inp_unit
    type(boltzmann_input_t) :: b

    !Namelist for boltzmann input
    NAMELIST /boltzmann_nml/ time_step,energy_eV,n_particles,boole_squared_moments,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, &
    & boole_antithetic_variate,boole_linear_temperature_simulation,i_integrator_type,seed_option

    open(newunit = b_inp_unit, file='boltzmann.inp', status='unknown')
    read(b_inp_unit,nml=boltzmann_nml)
    close(b_inp_unit)

    b%time_step = time_step
    b%energy_eV = energy_eV
    b%n_particles = n_particles
    b%density = density
    b%boole_squared_moments = boole_squared_moments
    b%boole_point_source = boole_point_source
    b%boole_collisions = boole_collisions
    b%boole_precalc_collisions = boole_precalc_collisions
    b%boole_refined_sqrt_g = boole_refined_sqrt_g
    b%boole_boltzmann_energies = boole_boltzmann_energies
    b%boole_linear_density_simulation = boole_linear_density_simulation
    b%boole_antithetic_variate = boole_antithetic_variate
    b%boole_linear_temperature_simulation = boole_linear_temperature_simulation
    b%i_integrator_type = i_integrator_type
    b%seed_option = seed_option
    b%num_particles = int(n_particles)

    print *,'GORILLA: Loaded input data from boltzmann.inp'

end subroutine read_boltzmann_inp_into_type

subroutine calc_boltzmann

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg,pi,echarge,ame,amp,clight
    use tetra_physics_mod, only: particle_mass,particle_charge,cm_over_e, coord_system, tetra_physics
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_set_num_threads
    use tetra_grid_settings_mod, only: n_field_periods, grid_size, grid_kind
    use tetra_grid_mod, only: ntetr, nvert, verts_rphiz, tetra_grid
    use gorilla_settings_mod, only: ispecies
    use collis_ions, only: stost
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use boltzmann_types_mod, only: filenames_t, output, moment_specs_t, counter_t, poloidal_flux_t, collisions_t, &
                                   boltzmann_input_t, iunits_t, time_t, start_t
    use boltzmann_writing_data_mod, only: write_data_to_files, name_files, unlink_files, open_files_before_main_loop, close_files
    use main_routine_supporting_functions_mod, only: print_progress, handle_lost_particles, initialise_loop_variables, &
    add_local_tetr_moments_to_output, normalise_prism_moments_and_prism_moments_squared, set_moment_specifications, &
    initialise_output, fourier_transform_moments, add_local_counter_to_counter, &
    get_ipert, calc_poloidal_flux, calc_starting_conditions, calc_collision_coefficients_for_all_tetrahedra, &
    carry_out_collisions

    integer :: kpart,i,j,n,l,m,k,p,ind_tetr,iface,iantithetic, i_part, n_prisms
    real(dp) :: v0,vpar,vperp, v
    real(dp), dimension(3) :: x
    real(dp), dimension(:,:), allocatable :: verts
    complex(dp), dimension(:,:), allocatable :: local_tetr_moments
    logical :: boole_initialized,boole_particle_lost
    type(start_t) :: start
    type(filenames_t) :: filenames
    type(moment_specs_t) :: moment_specs
    type(counter_t) :: counter, local_counter
    type(poloidal_flux_t) :: pflux
    type(collisions_t) :: c
    type(boltzmann_input_t) :: b
    type(iunits_t) :: iunits
    type(time_t) :: t

    call read_boltzmann_inp_into_type(b)
    call get_ipert()
    call initialize_gorilla(i_option,ipert)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete this again afterwards !!!!!!!!!!!!!!!!!!!!!!!
    if (ispecies.eq.4) particle_charge = 15*echarge
    print*, 'particle charge number = ', particle_charge/echarge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    v0=sqrt(2.d0*b%energy_eV*ev2erg/particle_mass)
    n_prisms = ntetr/3
    kpart = 0
    iantithetic = 1
    if (b%boole_antithetic_variate) iantithetic = 2

    call set_moment_specifications(moment_specs, b%boole_squared_moments)
    call initialise_output(moment_specs)
    call calc_square_root_g
    call calc_volume_integrals(b%boole_boltzmann_energies,b%boole_refined_sqrt_g, b%density, b%energy_eV)
    call calc_starting_conditions(b,v0,start,verts)
    call calc_poloidal_flux(pflux, verts)
    call calc_collision_coefficients_for_all_tetrahedra(b,c,v0,b%energy_eV)
    call name_files(filenames)
    call unlink_files(filenames)
    call open_files_before_main_loop(iunits)

    allocate(local_tetr_moments(moment_specs%n_moments,ntetr))
    if (b%i_integrator_type.eq.2) print*, 'Error: i_integratpr_type set to 2, this module only works with &
                                    & i_integrator_type set to 1'

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(counter,kpart,v0,start, iunits, moment_specs, pflux, b, c, iantithetic) &
    !$OMP& PRIVATE(p,l,n,boole_particle_lost,x,vpar,vperp,boole_initialized,ind_tetr,iface,t,v, local_tetr_moments,i,local_counter)
    print*, 'get number of threads', omp_get_num_threads()
    !$OMP DO

    !Loop over particles
    do p = 1,b%num_particles/iantithetic
        do l = 1,iantithetic

            n = (p-1)*iantithetic+l
            !$omp critical
            kpart = kpart+1 !in general not equal to n becuase of parallelisation
            call print_progress(b%num_particles,kpart,n)
            !$omp end critical

            call initialise_loop_variables(b, l, n, v0, start, local_counter,boole_particle_lost,t%step, t%confined, &
                                           local_tetr_moments,x,v,vpar,vperp)

            i = 0
            do while (t%confined.lt.b%time_step)
                i = i+1
                if (i.eq.1) then
                    boole_initialized = .false.
                endif

                if (b%boole_collisions) then
                    call carry_out_collisions(b, c, i, n, v0, t, x, vpar,vperp,ind_tetr, iface)
                    t%step = t%step/v0 !in carry_out_collisions, t%step is initiated as a length, so you need to divide by v0
                endif

                call orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t%step,boole_initialized,ind_tetr,iface,n, start,&
                            & local_tetr_moments, moment_specs, b, iunits, pflux,local_counter,t%remain)

                t%confined = t%confined + t%step - t%remain

                if (ind_tetr.eq.-1) then
                    call handle_lost_particles(iunits%et,t%confined, x, n, local_counter, boole_particle_lost)
                    exit
                endif

                v = sqrt(vpar**2+vperp**2)
            enddo

            !$omp critical
            counter%integration_steps = counter%integration_steps + i
            c%maxcol = max(dble(i)/dble(c%randcoli),c%maxcol)
            call add_local_counter_to_counter(local_counter,counter)
            !$omp end critical
            if (t%confined.eq.b%time_step) then
                write(iunits%rp,*) x, v, vpar, vperp, i, n
            endif
        enddo

        !$omp critical
        call add_local_tetr_moments_to_output(local_tetr_moments, moment_specs)
        !$omp end critical
    enddo !n
    !$OMP END DO
    !$OMP END PARALLEL

    call normalise_prism_moments_and_prism_moments_squared(moment_specs,b)
    if (moment_specs%n_moments.gt.0) call fourier_transform_moments(moment_specs)
    call close_files(iunits)
    call write_data_to_files(filenames,moment_specs)

    if (b%boole_precalc_collisions) print*, "maxcol = ", c%maxcol
    print*, 'Number of lost particles',counter%lost_particles
    print*, 'average number of pushings = ', counter%tetr_pushings/b%n_particles
    print*, 'average number of toroidal revolutions = ', counter%phi_0_mappings/b%n_particles
    print*, 'average number of integration steps = ', counter%integration_steps/b%n_particles
    PRINT*, 'particle mass = ', particle_mass
    PRINT*, 'absolute value of velocity = ', v0
    PRINT*, 'particle charge = ', particle_charge
    PRINT*, 'temperature = ', ev2erg*b%energy_eV
    print*, 'energy in eV = ', b%energy_eV
    print*, 'tracing time in seconds = ', b%time_step
    if((grid_kind.eq.2).or.(grid_kind.eq.3)) then
         print*, 'number of particles left through the outside = ', counter%lost_outside
         print*, 'number of particles left through the inside = ', counter%lost_inside
    endif

end subroutine calc_boltzmann

subroutine orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, n,start,local_tetr_moments, &
                                            moment_specs, b, iunits, pflux, local_counter, t_remain_out)

    use pusher_tetra_rk_mod, only: pusher_tetra_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type, boole_array_optional_quantities
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: vperp_func
    use find_tetra_mod, only: find_tetra
    use tetra_grid_mod, only: tetra_grid, ntetr
    use boltzmann_types_mod, only: moment_specs_t, counter_t, poloidal_flux_t, boltzmann_input_t, iunits_t, start_t
    use tetra_grid_settings_mod, only: grid_kind
    use orbit_timestep_gorilla_supporting_functions_mod, only: calc_and_write_poincare_mappings_and_divertor_intersections, &
            categorize_lost_particles, update_local_tetr_moments, initialize_constants_of_motion, calc_particle_weights_and_jperp

    type(moment_specs_t), intent(in)             :: moment_specs
    type(iunits_t), intent(in)                   :: iunits
    type(counter_t), intent(inout)               :: local_counter
    type(poloidal_flux_t)                        :: pflux
    type(boltzmann_input_t)                      :: b
    type(start_t)                                :: start
    real(dp), dimension(3), intent(inout)        :: x
    complex(dp), dimension(:,:), intent (inout)  :: local_tetr_moments
    real(dp), intent(inout)                      :: vpar,vperp
    real(dp), intent(in)                         :: t_step
    logical, intent(inout)                       :: boole_initialized
    integer, intent(inout)                       :: ind_tetr,iface
    real(dp), intent(out), optional              :: t_remain_out
    real(dp), dimension(3)                       :: z_save, x_save
    real(dp)                                     :: t_remain,t_pass,perpinv, aphi
    logical                                      :: boole_t_finished
    integer                                      :: ind_tetr_save,iper_phi,n
    type(optional_quantities_type)               :: optional_quantities
    
    if(.not.boole_initialized) then !If orbit_timestep is called for the first time without grid position
        call check_coordinate_domain(x) !Check coordinate domain (optionally perform modulo operation)
        call find_tetra(x,vpar,vperp,ind_tetr,iface) !Find tetrahedron index and face index for position x
        if(ind_tetr.eq.-1) then !If particle doesn't lie inside any tetrahedron
            t_remain_out = t_step
            return
        endif
        z_save = x-tetra_physics(ind_tetr)%x1
        call calc_particle_weights_and_jperp(b,n,z_save,vpar,vperp,ind_tetr,pflux,start)
        boole_initialized = .true.
    endif
          
    if(t_step.eq.0.d0) return !Exit the subroutine after initialization, if time step equals zero
    if(boole_initialized) z_save = x-tetra_physics(ind_tetr)%x1
    call initialize_constants_of_motion(vperp,z_save,ind_tetr,perpinv)
    t_remain = t_step
    boole_t_finished = .false.
    local_counter%tetr_pushings = local_counter%tetr_pushings -1 !set tetr_pushings to -1 because when entering the loop it will go back to one without pushing

    do !Loop for tetrahedron pushings until t_step is reached
        local_counter%tetr_pushings = local_counter%tetr_pushings +1

        if(ind_tetr.eq.-1) then
            if((grid_kind.eq.2).or.(grid_kind.eq.3)) then
                call categorize_lost_particles(ind_tetr_save,z_save,pflux,local_counter,t_remain,t_remain_out)
            endif
            exit
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

        call calc_and_write_poincare_mappings_and_divertor_intersections(x_save,x,n,iper_phi,local_counter,iunits,boole_t_finished)
        call update_local_tetr_moments(moment_specs,local_tetr_moments,ind_tetr_save,n,start,optional_quantities)
        if(boole_t_finished) then !Orbit stops within cell, because "flight"-time t_step has finished
            if( present(t_remain_out)) t_remain_out = t_remain
            exit
        endif

    enddo !Loop for tetrahedron pushings

    vperp = vperp_func(z_save,perpinv,ind_tetr_save) !Compute vperp from position

end subroutine orbit_timestep_gorilla_boltzmann

end module boltzmann_mod