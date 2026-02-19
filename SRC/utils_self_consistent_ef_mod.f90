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

subroutine parallelised_particle_pushing(species,j,boole_diffusion_coefficient,n_particles_in)

    use gorilla_applets_types_mod, only: counter, c, in, time_t, moment_specs, counter_t, particle_status_t, start, exit_data, s, &
    g, maximum_s
    use tetra_grid_mod, only: ntetr
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use utils_parallelised_particle_pushing_mod, only: print_progress, handle_lost_particles, add_local_tetr_moments_to_output, &
    add_local_counter_to_counter, initialise_loop_variables, carry_out_collisions, update_exit_data, update_start_type, &
    initialise_seed_for_random_numbers_for_each_thread
    use russian_roulette_mod, only: play_russian_roulette, rr, local_rr_t, initiate_local_rr
    use constants, only: ev2erg
    

    integer, intent(in)                               :: species, j
    integer, intent(in), optional                     :: n_particles_in
    logical, intent(in)                               :: boole_diffusion_coefficient
    integer                                           :: kpart, iantithetic, ind_tetr, iface, n_particles
    integer                                           :: p, l, n, i, k
    real(dp), dimension(3)                            :: x
    real(dp)                                          :: vpar,vperp, t_step_s, v, v_save,vpar_save, vperp_save, t_tot, v_init
    type(time_t)                                      :: t
    type(counter_t)                                   :: local_counter
    type(particle_status_t)                           :: particle_status
    complex(dp), dimension(:,:), allocatable          :: local_tetr_moments
    logical                                           :: thread_flag = .true.
    type(local_rr_t)                                  :: local_rr
    real(dp), dimension(10)                           :: particle_state_for_rr

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

    maximum_s = 0.0_dp
    

    allocate(local_tetr_moments(moment_specs%n_moments,ntetr))
    kpart = 0
    iantithetic = 1
    if (in%boole_antithetic_variate) iantithetic = 2

    if (boole_diffusion_coefficient) then
        s%delta_s = 0.0_dp
        s%delta_s_squared = 0.0_dp
        s%check = 0
        s%f_v = 0
    endif

    t_tot = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(counter, kpart,species, in, c, iantithetic, start, j, s, boole_diffusion_coefficient,n_particles,rr) &
    !$OMP& REDUCTION(+:t_tot) &
    !$OMP& PRIVATE(p,l,n,i,x,v_save,vpar,vperp,t,ind_tetr,iface,local_tetr_moments,local_counter,particle_status,t_step_s,k,v, &
    !$OMP& vpar_save, vperp_save, particle_state_for_rr, v_init) &
    !$OMP& FIRSTPRIVATE(thread_flag,local_rr)
    if (omp_get_thread_num().eq.0) print*, 'get number of threads', omp_get_num_threads()
    !$OMP DO SCHEDULE(static)
    !SCHEDULE(dynamic,1)
    !SCHEDULE(static)

    !Loop over particles
    do p = 1,n_particles/iantithetic

        if ((.not.in%boole_precalc_collisions).and.thread_flag) then
            call initialise_seed_for_random_numbers_for_each_thread(omp_get_thread_num(), j)
            thread_flag = .false.
        endif

        do l = 1,iantithetic
            n = (p-1)*iantithetic+l
            !$omp atomic update
            kpart = kpart + 1 !in general not equal to n because of parallelisation
            call print_progress(n_particles,kpart,n)

            call initialise_loop_variables(l, n, local_counter,particle_status,t,local_tetr_moments,x,vpar,vperp,species)

            i = 0

            if (rr%boole_russian_roulette) call initiate_local_rr(local_rr,100)

    ! if (s%s0.lt.2.0d-2) then
    !      open(1000+n)
    ! endif
            v_init = sqrt(vpar**2+vperp**2)

            do while (t%confined.lt.start%t(species))

                i = i+1

                vpar_save = vpar
                vperp_save = vperp
                v_save = sqrt(vpar_save**2+vperp_save**2)

                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, t, x, vpar,vperp,ind_tetr, iface, species)
                    t%step = t%step/start%v0(species) !in carry_out_collisions, t%step is initiated as a length, so you need to divide by v0
                endif

                v = sqrt(vpar**2+vperp**2)

                if (rr%boole_russian_roulette.and.(i.gt.1)) then
                    if (local_rr%boole_eliminated) then
                        call initiate_next_split_particle(local_rr,vpar,vperp,t,x,ind_tetr,iface,particle_status,n,species)
                        if (local_rr%boole_eliminated) exit
                    else
                        particle_state_for_rr = (/vpar,vperp,t%confined,t%remain,t%step,x,dble(ind_tetr),dble(iface)/)
                        call play_russian_roulette(start%weight(n,species),v,v_save,particle_state_for_rr,local_rr)
                    endif
                endif

                if (boole_diffusion_coefficient) then
                    t_step_s = start%t(species)/s%k - (t%confined - start%t(species)/s%k*int(t%confined/(start%t(species)/s%k)))
                    k = int(t%confined/(start%t(species)/s%k)) + 1
                else 
                    k = 1.0_dp !if boole_diffusion_coefficient.eqv..false., k is meaningless
                endif

                call orbit_timestep_gorilla_self_consistent_ef(x,vpar,vperp,t,particle_status,ind_tetr,iface,n,&
                &local_tetr_moments, local_counter, species, j, t_step_s, k, boole_diffusion_coefficient,local_rr)

                t%confined = t%confined + t%step - t%remain
                t_tot = t_tot + t%step - t%remain



                if ((ind_tetr.eq.-1).and.(.not.rr%boole_russian_roulette)) then
                    call handle_lost_particles(local_counter, particle_status%lost)
                    exit
                endif
            enddo

            !print*, 'Confinement time for particle ', n, ' (species ', species, '): ', t%confined
            !print*, 'energy change for particle', n, ' (species ', species, '):', &
            !v**2*start%particle_mass(species)/(ev2erg*2*start%energy(n,species))

    ! if (s%s0.lt.2.0d-2) then
    !      close(1000+n)
    ! endif

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

    print*, 'Total tracing time of all particles divided by number of particles is: ', t_tot/n_particles, 's'

end subroutine parallelised_particle_pushing

subroutine orbit_timestep_gorilla_self_consistent_ef(x,vpar,vperp,t,particle_status,ind_tetr,iface, n,local_tetr_moments, &
                                local_counter,species, j, t_step_s, k, boole_diffusion_coefficient,local_rr)

    use pusher_tetra_rk_mod, only: pusher_tetra_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: vperp_func
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_types_mod, only: counter_t, particle_status_t, g, start, s, in, time_t, maximum_s
    use tetra_grid_settings_mod, only: grid_kind, sfc_s_min
    use utils_orbit_timestep_mod, only: identify_particles_entering_annulus, update_local_tetr_moments, &
                                        initialize_constants_of_motion, compute_radial_fluxes
    use constants, only: echarge
    use russian_roulette_mod, only: play_russian_roulette, rr, local_rr_t

    integer, intent(in)                          :: species, j, n
    logical, intent(in)                          :: boole_diffusion_coefficient
    real(dp), intent(inout)                      :: t_step_s
    type(counter_t), intent(inout)               :: local_counter
    type(particle_status_t), intent(inout)       :: particle_status
    type(time_t)                                 :: t
    real(dp), dimension(3), intent(inout)        :: x
    complex(dp), dimension(:,:), intent (inout)  :: local_tetr_moments
    real(dp), intent(inout)                      :: vpar,vperp
    integer, intent(inout)                       :: ind_tetr,iface
    real(dp), dimension(3)                       :: z_save, x_save, z_save_at_x_save
    real(dp)                                     :: t_pass,perpinv, vpar_save, vperp_save
    logical                                      :: boole_t_finished, boole_lost_inside, mirror_condition
    integer                                      :: ind_tetr_save,iper_phi,k,i
    type(optional_quantities_type)               :: optional_quantities
    real(dp)                                     :: tau, v, t_pusher, v_save, vpar_init, vperp_init, critical_distance
    type(local_rr_t)                             :: local_rr
    real(dp), dimension(10)                      :: particle_state_for_rr


    v = sqrt(vpar**2+vperp**2)
    if(.not.particle_status%initialized) then !If orbit_timestep is called for the first time without grid position
        call check_coordinate_domain(x) !Check coordinate domain (optionally perform modulo operation)
        call find_tetra(x,vpar,vperp,ind_tetr,iface) !Find tetrahedron index and face index for position x
        if(ind_tetr.eq.-1) then !If particle doesn't lie inside any tetrahedron
            t%remain = t%step
            return
        endif
        z_save = x-tetra_physics(ind_tetr)%x1
        if ((j.eq.1).or.(.not.in%boole_static_ne)) then
            call calc_particle_weights_and_jperp(n,z_save,vpar,vperp,ind_tetr,species,boole_diffusion_coefficient)
        endif
        if ((.not.in%boole_static_ne).and.(species.eq.1)) then
            start%x(:,n,2) = start%x(:,n,1)
            start%weight(n,2) = start%weight(n,1)
        endif 
        particle_status%initialized = .true.
    endif
          
    if(t%step.eq.0.0_dp) return !Exit the subroutine after initialization, if time step equals zero
    if(particle_status%initialized) z_save = x-tetra_physics(ind_tetr)%x1
    call initialize_constants_of_motion(vperp,z_save,ind_tetr,perpinv)
    t%remain = t%step
    boole_t_finished = .false.
    local_counter%tetr_pushings = local_counter%tetr_pushings -1 !set tetr_pushings to -1 because when entering the loop it will go back to one without pushing
    vpar_init = vpar
    vperp_init = vperp

    do !Loop for tetrahedron pushings until t%step is reached
        local_counter%tetr_pushings = local_counter%tetr_pushings +1

        if(ind_tetr.eq.-1) then
            mirror_condition = .not.(x(1).gt.1.01_dp*sfc_s_min)
            if (in%boole_static_ne) mirror_condition = (.not.(x(1).gt.1.01_dp*sfc_s_min)).or.(.not.(x(1).lt.0.99_dp))
            if(mirror_condition) then !.or.(.not.(x(1).lt.0.99_dp)) <-- include this if you also want a mirror term at s=1
                call mirror_particles_on_domain_boundaries(x,vpar,n,ind_tetr,iface,z_save,perpinv,ind_tetr_save)
                if (ind_tetr.eq.-1) exit
            elseif (x(1).lt.0.99_dp) then
                call treat_particles_that_are_lost_but_should_not_be(z_save_at_x_save,ind_tetr_save,z_save,x_save,x,vpar,vperp, &
                                                                     perpinv,ind_tetr,vpar_save,vperp_save,vpar_init,vperp_init)
            else
                exit
            endif
        endif

        ind_tetr_save = ind_tetr
        x_save = x
        z_save_at_x_save = z_save
        vpar_save = vpar
        vperp_save = vperp!_func(z_save_at_x_save,perpinv,ind_tetr_save)
        v_save = sqrt(vpar_save**2+vperp**2)
        t_pusher = t%remain
        if (boole_diffusion_coefficient) t_pusher = min(t%remain,t_step_s)

        select case(ipusher) !Calculate trajectory
            case(1)
                call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t%remain,t_pass,boole_t_finished,iper_phi)
            case(2)
                call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_pusher, &
                                                    & t_pass,boole_t_finished,iper_phi,optional_quantities)
        end select

        vperp = vperp_func(z_save,perpinv,ind_tetr_save) !Compute vperp from position
        v = sqrt(vperp**2+vpar**2)

        if (boole_diffusion_coefficient) then
            if (boole_t_finished.and.(t%remain.ge.t_step_s)) then
                if (t%remain.gt.t_step_s) boole_t_finished = .false.
                t_step_s = start%t(species)/s%k + t_pass
                local_counter%tetr_pushings = local_counter%tetr_pushings -1
                critical_distance = 1.0_dp-s%s0
                if (abs(x(1)-s%s0).gt.critical_distance) s%boole_large_distance(n)=.true.
                
                !$omp critical
                s%delta_s(k) = s%delta_s(k) + (x(1) - s%s0)*start%weight(n,species)
                s%delta_s_squared(k) = s%delta_s_squared(k) + (x(1) - s%s0)**2*start%weight(n,species)
                if (s%boole_large_distance(n)) then !delete contribution to delta_s, double contribution to delta_s_squared
                    s%delta_s(k) = s%delta_s(k) - (x(1) - s%s0)*start%weight(n,species)
                    s%delta_s_squared(k) = s%delta_s_squared(k) + (x(1) - s%s0)**2*start%weight(n,species)
                endif
                s%check(k) = s%check(k) + 1
                i = min(int(s%j/10*v/start%v0(species))+1, s%j)
                if (int(10*v/start%v0(species))+1.gt.s%j) print*, 'ATTENTION: particle is faster than 10*v_t'
                s%f_v(k,i) = s%f_v(k,i) + (x(1) - s%s0)**2/s%n_particles 
                !if (n.eq.1520) write(71,*) s%time(k), x, vpar, vperp
                !$omp end critical
                ! if (s%s0.lt.2.0d-2) then
                !      write(1000+n,*) s%time(k), x, vpar, v
                ! endif
                k = k+1
            endif
        endif

        t%remain = t%remain - t_pass
        t_step_s = t_step_s - t_pass

        call update_local_tetr_moments(local_tetr_moments,ind_tetr_save,n,optional_quantities,species)
        if((grid_kind.eq.2).or.(grid_kind.eq.3)) call compute_radial_fluxes(ind_tetr_save,ind_tetr,x)

        if (rr%boole_russian_roulette) then
            particle_state_for_rr = (/vpar,vperp,t%confined,t%remain,t%step,x,dble(ind_tetr),dble(iface)/)
            call play_russian_roulette(start%weight(n,species),v,v_save,particle_state_for_rr,local_rr)
        endif

        if(boole_t_finished.or.local_rr%boole_eliminated) then !Orbit stops within cell, because "flight"-time t%step has finished
            exit
        endif

        ! !$omp critical
        ! if (x(1).gt.maximum_s) then
        !     !print*, 'maximum s-value increased to ', x(1)
        !     maximum_s = x(1)
        ! endif
        ! !$omp end critical

    enddo !Loop for tetrahedron pushings

end subroutine orbit_timestep_gorilla_self_consistent_ef

subroutine initiate_next_split_particle(local_rr,vpar,vperp,t,x,ind_tetr,iface,particle_status,n,species)

    use gorilla_applets_types_mod, only : particle_status_t, time_t, start
    use russian_roulette_mod, only: local_rr_t, prepare_next_split_particle

    type(local_rr_t):: local_rr
    real(dp) :: vpar,vperp
    integer :: ind_tetr, iface, n, species
    type(time_t) :: t
    real(dp), dimension(3) :: x
    type(particle_status_t) :: particle_status
    integer :: id

    call prepare_next_split_particle(local_rr,id)

    if (local_rr%boole_eliminated.eqv..false.) then

        vpar =                    local_rr%particle_state(id,1)
        vperp =                   local_rr%particle_state(id,2)
        t%confined =              local_rr%particle_state(id,3)
        t%remain =                local_rr%particle_state(id,4)
        t%step =                  local_rr%particle_state(id,5)
        x =                       local_rr%particle_state(id,6:8)
        ind_tetr =                int(local_rr%particle_state(id,9))
        iface =                   int(local_rr%particle_state(id,10))
        start%weight(n,species) = local_rr%weight(id)

        particle_status%lost = .false.
        particle_status%initialized = .true.
        particle_status%exit = .false.
    endif

end subroutine initiate_next_split_particle

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

    !boole_keep_phi_elec=.true. prevents make_tetra_physics from overwriting phi_elec
    call make_tetra_physics(coord_system,ipert,boole_keep_phi_elec=.true.)

    if (i.gt.0) call print_data(i)

end subroutine perform_electric_potential_update

subroutine calc_phi_elec_from_rho(i)

    use gorilla_applets_types_mod, only: in, ep, start, output
    use constants, only: ev2erg, eps, echarge

    integer, intent(in) :: i
    real(dp) :: factor, factor_from_tracing_time

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

        !decrease factor in case very few particles are simulated
        !if (in%n_particles.lt.100.0_dp) factor = factor*in%n_particles/100.0_dp

        
        ep%phi_elec_from_rho =  ep%rho_vert*factor

        factor_from_tracing_time = 10.0_dp
        !tracing time is about 25 transport times, since densities are normalised by tracing time but particles only stay inside the 
        !computation domain for about the transport time, this factor is compensated here (with some safety margin)
        if (.not.in%boole_static_ne) ep%phi_elec_from_rho = ep%phi_elec_from_rho*factor_from_tracing_time
        
        !if (i.gt.0) ep%phi_elec_from_rho = ep%phi_elec_from_rho/sqrt(dble(i))

        ! ep%phi_elec_from_rho = ep%phi_elec_from_rho*(3.5d3*1.6022d-12/4.8032d-10)/maxval(abs(ep%phi_elec_from_rho))
        ! if (i.gt.1) ep%phi_elec_from_rho=0
        ! print*, 'maximum phi is', maxval(ep%phi_elec_from_rho), 'minimum phi is', minval(ep%phi_elec_from_rho)
        
    endif

end subroutine calc_phi_elec_from_rho

subroutine calc_average_charge_density_per_flux_layer(i)

    use gorilla_applets_types_mod, only:  g, output, ep, in, start, one_d, exit_data
    use tetra_grid_mod, only: nvert, verts_sthetaphi
    use tetra_grid_settings_mod, only: grid_size
    use constants, only: echarge

    integer, intent(in) :: i 
    integer :: ns, j, species, ind_prism, n_species
    real(dp), dimension(:), allocatable :: electron_densities
    real(dp) :: s, electron_density_factor, factor_from_ion_weights

    n_species = in%n_species

    allocate(electron_densities(grid_size(1)))
    if (one_d%boole_print_densities) then
        if (.not.allocated(one_d%densities)) allocate(one_d%densities(grid_size(1),in%n_species))
        one_d%densities = 0.0_dp
    endif

    do j = 1,in%num_particles
        do species = 1,2
            if (i.gt.0) ep%total_tracing_time(i) = ep%total_tracing_time(i) + exit_data%t_confined(j,species)
        enddo
    enddo

    !If in%boole_static_ne compute electron_densities
    if (in%boole_static_ne) then
        do ns = 1,grid_size(1)
            s = (verts_sthetaphi(1, grid_size(3)*(ns-1)+1) + verts_sthetaphi(1,grid_size(3)*ns+1))/2.0_dp
            electron_densities(ns) = in%density*(1.0_dp-0.9_dp*s)
        enddo
    endif

    if (in%boole_static_ne) then
        n_species = in%n_species-1
        if (in%boole_linear_density_simulation) then
            electron_density_factor = 1.0d0
        else
            factor_from_ion_weights = sum(start%weight(:,1))/(in%num_particles*in%density*sum(output%prism_volumes(:)))
            electron_density_factor = in%density/(sum(electron_densities*ep%s_shell_volumes)/sum(ep%s_shell_volumes))*&
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
        ep%rho_flux_layer(ns) = ep%rho_flux_layer(ns)/ep%s_shell_volumes(ns)
        if (one_d%boole_print_densities) one_d%densities(ns,:) = one_d%densities(ns,:)/ep%s_shell_volumes(ns)
    enddo

end subroutine calc_average_charge_density_per_flux_layer

subroutine calc_rho_on_vertices

    use tetra_grid_mod, only: nvert, verts_rphiz, verts_sthetaphi
    use tetra_grid_settings_mod, only: grid_size
    use gorilla_applets_types_mod, only: g, ep, in
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
        !value_to_be_set = -sin(dble(ns-1)/dble(grid_size(1))*pi)
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
    if (.not.in%boole_static_ne) value_to_be_set = -sum(ep%rho_flux_layer*ep%s_shell_volumes)/sum(ep%s_shell_volumes)
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
                                                            perpinv,ind_tetr,vpar_save,vperp_save,vpar_init,vperp_init)

    use utils_orbit_timestep_mod, only: initialize_constants_of_motion
    use supporting_functions_mod, only: vperp_func

    integer, intent(in)                          :: ind_tetr_save
    real(dp), dimension(3), intent(in)           :: z_save_at_x_save, x_save
    real(dp), dimension(3), intent(inout)        :: x, z_save
    real(dp), intent(inout)                      :: vpar,vperp, perpinv, vpar_save, vperp_save, vpar_init,vperp_init
    integer, intent(inout)                       :: ind_tetr
    real(dp)                                     :: v
    integer                                      :: problem_unit

    print*, 'This should not happen.'
    vperp      = vperp_func(z_save,          perpinv,ind_tetr_save)

    print*, 'x_save, vpar_save, vperp_save = ', x_save, vpar_save, vperp_save
    print*, 'x, vpar, vperp = ', x, vpar, vperp
    print*, 'vpar/vpar_init, vperp/vperp_init = ', vpar/vpar_init, vperp/vperp_init

    !$omp critical
    open(newunit = problem_unit, file = 'pushing_problems.dat', position = 'append')
    write(problem_unit,*) 'x_save, vpar_save, vperp_save = ', x_save, vpar_save, vperp_save
    write(problem_unit,*) 'x, vpar, vperp = ', x, vpar, vperp
    close(problem_unit)
    !$omp end critical

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

subroutine calc_particle_weights_and_jperp(n,z_save,vpar,vperp,ind_tetr, species, boole_diffusion_coefficient)

    use gorilla_applets_types_mod, only: in, start, s
    use tetra_physics_mod, only: tetra_physics
    use constants, only: ev2erg, pi
    use volume_integrals_and_sqrt_g_mod, only: sqrt_g
    use supporting_functions_mod, only: bmod_func

    logical, intent(in) :: boole_diffusion_coefficient
    real(dp), intent(in) :: vpar, vperp
    real(dp), dimension(3), intent(in) :: z_save
    real(dp), dimension(3) :: x
    integer, intent(in) :: n,ind_tetr
    integer, intent(in) :: species
    real(dp) :: phi_elec_func, temperature

    x = tetra_physics(ind_tetr)%x1 + z_save
    start%weight(n,species) = start%weight(n,species)*abs((tetra_physics(ind_tetr)%sqg1 + sum(tetra_physics(ind_tetr)%gsqg*z_save)))
    !print*, 'weight before = ', start%weight(n,species), n

    if (in%boole_linear_density_simulation) then
        start%weight(n,species) = start%weight(n,species)*(1.0_dp-0.9_dp*x(1))
    endif

    if (boole_diffusion_coefficient) then
        start%weight(n,species) = 1.0_dp/s%n_particles
    elseif (in%boole_boltzmann_energies) then

        temperature = in%energy_eV
        if (boole_diffusion_coefficient) temperature = s%temperature
        start%weight(n,species) = start%weight(n,species)/  &
        (sqrt(start%energy(n,species)/temperature)*exp(-start%energy(n,species)/temperature)/(temperature*ev2erg*sqrt(pi)/2))
        !the last term is the integral of the function from zero to inf over energy in correct units
        !start%weight(n,species) = start%weight(n,species)/((1.0_dp+(start%energy(n,species)/s%temperature)**(3.5_dp))*  &
        !sqrt(start%energy(n,species)/s%temperature)*exp(-start%energy(n,species)/s%temperature)/(s%temperature*ev2erg*14035.0_dp))
        !the last term is the integral of the function from zero to inf over energy in correct units

        phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z_save)
        phi_elec_func = 0.0_dp !when working with fixed sources, electric potentials are not useful since they change the weights
        !and thus the magnitude o the sources

        start%weight(n,species) = start%weight(n,species)*2/sqrt(pi)*sqrt(start%energy(n,species)*ev2erg)
        ! if (.not.boole_diffusion_coefficient) then 
        !     start%weight(n,species) = start%weight(n,species)*start%epsilon_max*in%energy_eV*ev2erg
        ! endif
        if (.not. in%boole_linear_temperature_simulation) then
            temperature = in%energy_eV*ev2erg
            if (boole_diffusion_coefficient) temperature = s%temperature*ev2erg
            start%weight(n,species) =start%weight(n,species)/temperature**1.5_dp* &
                        & exp(-(start%energy(n,species)*ev2erg+start%particle_charge(species)*phi_elec_func)/temperature)
        else
            temperature = in%energy_eV*ev2erg*(1.0_dp-0.9*x(1))
            if (boole_diffusion_coefficient) temperature = s%temperature*ev2erg*(1.0_dp-0.9*x(1))
            start%weight(n,species) = start%weight(n,species)/temperature**1.5_dp* &
                        & exp(-(start%energy(n,species)*ev2erg+start%particle_charge(species)*phi_elec_func)/temperature)
        endif
    endif

    start%jperp(n,species) = start%particle_mass(species)*vperp**2*start%cm_over_e(species)/(2*bmod_func(z_save,ind_tetr))*(-1)

    
    !-1 because of negative gyrophase
    
!print*, 'weight after = ', start%weight(n,species)*(1.0_dp+(start%energy(n,species)/s%temperature)**(3.5_dp))/14035.0_dp/2*sqrt(pi)
!print*, 'weight after for real = ', start%weight(n,species)



end subroutine calc_particle_weights_and_jperp

subroutine calc_starting_conditions

    use gorilla_applets_types_mod, only: in
    use tetra_grid_settings_mod, only: sfc_s_min
    
    real(dp), dimension(:,:,:), allocatable                :: rand_matrix

    allocate(rand_matrix(5,in%num_particles,in%n_species))
    call RANDOM_NUMBER(rand_matrix)

    call allocate_start_type
    call set_particle_type_specifications
    call set_starting_positions(rand_matrix)!,s0=sfc_s_min*1.1_dp)
    call set_rest_of_individual_particle_specifications(rand_matrix)

end subroutine calc_starting_conditions

subroutine allocate_start_type(n_particles_in)

    use gorilla_applets_types_mod, only: start, in

    integer, intent(in), optional :: n_particles_in
    integer :: n_particles

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

    !before allocating, deallocate if necessary
    call deallocate_start_type

    allocate(start%x(3,n_particles,in%n_species))
    allocate(start%pitch(n_particles,in%n_species))
    allocate(start%energy(n_particles,in%n_species))
    allocate(start%weight(n_particles,in%n_species))
    allocate(start%jperp(n_particles,in%n_species))
    allocate(start%lost(n_particles,in%n_species))
    allocate(start%particle_charge(in%n_species))
    allocate(start%particle_mass(in%n_species))
    allocate(start%cm_over_e(in%n_species))
    allocate(start%t(in%n_species))
    allocate(start%v0(in%n_species))

end subroutine allocate_start_type

subroutine deallocate_start_type

    use gorilla_applets_types_mod, only: start

    if (allocated(start%x))               deallocate(start%x)
    if (allocated(start%pitch))           deallocate(start%pitch)
    if (allocated(start%energy))          deallocate(start%energy)
    if (allocated(start%weight))          deallocate(start%weight)
    if (allocated(start%jperp))           deallocate(start%jperp)
    if (allocated(start%lost))            deallocate(start%lost)
    if (allocated(start%particle_charge)) deallocate(start%particle_charge)
    if (allocated(start%particle_mass))   deallocate(start%particle_mass)
    if (allocated(start%cm_over_e))       deallocate(start%cm_over_e)
    if (allocated(start%t))               deallocate(start%t)
    if (allocated(start%v0))              deallocate(start%v0)

end subroutine deallocate_start_type

subroutine set_starting_positions(rand_matrix,species_in,s0)

    use gorilla_applets_types_mod, only: in, start, g
    use tetra_physics_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind
    use constants, only: pi
    use tetra_grid_settings_mod, only: sfc_s_min, n_field_periods

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix
    integer, dimension(:), intent(in), optional :: species_in
    integer, dimension(:), allocatable :: species
    integer :: i
    real(dp), intent(in), optional :: s0

    if (present(species_in)) then 
        allocate(species(size(species_in)))
        species = species_in
    else
        allocate(species(in%n_species))
        species = [(i,i=1,in%n_species)]
    endif

    start%x(1,:,species) = sfc_s_min + rand_matrix(1,:,:)*(1-sfc_s_min)
    if (present(s0).and.(.not.in%boole_static_ne))  start%x(1,:,species) = s0
    start%x(2,:,species) = 2*pi*rand_matrix(2,:,:) !theta
    start%x(3,:,species) = 2*pi/n_field_periods*rand_matrix(3,:,:) !phi

    !unless a single species is initiated, make elctrons and ions start at identical positions in real space
    if (size(species).gt.1) then 
        start%x(:,:,in%n_species) = start%x(:,:,1)
    endif

end subroutine set_starting_positions

subroutine set_rest_of_individual_particle_specifications(rand_matrix,boole_diffusion_coefficient_in,species_in,n_particles_in)

    use gorilla_applets_types_mod, only: in, start, s

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix
    integer, dimension(:), intent(in), optional :: species_in
    logical, intent(in), optional :: boole_diffusion_coefficient_in
    logical :: boole_diffusion_coefficient=.false.
    integer, dimension(:), allocatable :: species
    integer :: i
    integer, intent(in), optional :: n_particles_in
    integer :: n_particles
    real(dp), dimension(:,:), allocatable :: radial_transport_energies
    real(dp) :: temperature

    if (present(species_in)) then 
        allocate(species(size(species_in)))
        species = species_in
    else
        allocate(species(in%n_species))
        species = [(i,i=1,in%n_species)]
    endif

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

    if (present(boole_diffusion_coefficient_in)) boole_diffusion_coefficient = boole_diffusion_coefficient_in

    start%pitch(:,species) = 2*rand_matrix(4,:,:)-1 !pitch parameter
    start%energy(:,species) = in%energy_eV
    if (in%boole_boltzmann_energies.or.boole_diffusion_coefficient) then
        !start%energy(:,species) = start%epsilon_max*in%energy_eV*rand_matrix(5,:,:) !boltzmann energy distribution
        temperature = in%energy_eV
        if (boole_diffusion_coefficient) temperature = s%temperature
        allocate(radial_transport_energies(n_particles,size(species)))
        call generate_distribution_sqrt_x_exp_neg_x(start%epsilon_max,radial_transport_energies)
        !call generate_marker_distribution(start%epsilon_max,radial_transport_energies)
        start%energy(:,species) = temperature*radial_transport_energies !boltzmann energy distribution
    endif
    
    if (in%boole_antithetic_variate) then
        start%x(:,1:n_particles:2,species) = start%x(:,2:n_particles:2,species)
        start%pitch(1:n_particles:2,species) = -start%pitch(2:n_particles:2,species)
        start%energy(1:n_particles:2,species) = start%energy(2:n_particles:2,species)
    endif

end subroutine set_rest_of_individual_particle_specifications

subroutine set_particle_type_specifications

    use gorilla_applets_types_mod, only: in, start
    use constants, only: echarge,ame,clight, ev2erg, amp
    use gorilla_settings_mod, only: ispecies

    real(dp) :: charge, mass, cm_over_e

    select case(ispecies)
        case(1) !electron
            charge = -echarge
            mass = ame
            cm_over_e = -clight*ame/echarge
        case(2) !deuterium ion
            charge = echarge
            mass = 2.d0*amp
            cm_over_e=2.d0*clight*amp/echarge
        case(3) !alpha particle
            charge = 2.d0*echarge
            mass = 4.d0*amp
            cm_over_e=2.d0*clight*amp/echarge
        case(4) !ionised tungsten
            charge = 74.d0*echarge
            mass = 184.d0*amp
            cm_over_e= 184.d0*clight*amp/(74.d0*echarge)
    end select  

    start%particle_charge = (/charge, -echarge/)
    start%particle_mass = (/mass, ame/)
    start%cm_over_e = (/cm_over_e, -clight*ame/echarge/)
    start%t = (/in%time_step, in%time_step/) !/42.0_dp
    if (in%boole_static_ne) start%t(2) = 0.0_dp

    start%v0 = sqrt(2.0_dp*in%energy_eV*ev2erg/start%particle_mass)
    start%epsilon_max = 16.0_dp

    call set_weight

end subroutine set_particle_type_specifications

subroutine set_weight

    use gorilla_applets_types_mod, only: start, in
    use tetra_grid_settings_mod, only: sfc_s_min, n_field_periods
    use constants, only: pi

    start%weight = in%density*(1-sfc_s_min)*4*pi**2/n_field_periods

end subroutine set_weight

subroutine calc_s_shell_volumes

    use gorilla_applets_types_mod, only:  g, output, ep
    use tetra_grid_mod, only: verts_sthetaphi
    use tetra_grid_settings_mod, only: grid_size

    integer :: ns, j, ind_prism
    real(dp) :: s

    allocate(ep%s_shell_volumes(grid_size(1)))
    ep%s_shell_volumes = 0.0_dp

    do ns = 1,grid_size(1)
        s = (verts_sthetaphi(1, grid_size(3)*(ns-1)+1) + verts_sthetaphi(1,grid_size(3)*ns+1))/2.0_dp
        do j = 1, grid_size(2)*grid_size(3)*2
            ind_prism = g%prisms_per_flux_tube(ns,j)
            ep%s_shell_volumes(ns) = ep%s_shell_volumes(ns) + output%prism_volumes(ind_prism)
        enddo
    enddo

end subroutine calc_s_shell_volumes

subroutine calc_electron_diffusion_coefficients !call this before the first ion pushing

    use gorilla_applets_types_mod, only: in, dc, start, s, g, exit_data
    use tetra_grid_settings_mod, only: grid_size
    use tetra_grid_mod, only: verts_sthetaphi
    use utils_data_pre_and_post_processing_mod, only: prepare_next_round_of_parallelised_particle_pushing, &
    calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared, initialize_exit_data
    use llsq_mod, only: llsq
    use tetra_physics_mod, only: particle_mass,tetra_physics
    use russian_roulette_mod, only: prepare_russian_roulette
    use tetra_grid_settings_mod, only: sfc_s_min

    integer :: ns, i, n_particles, file_id, n_ignored, ns_start
    real(dp) :: extrapolation_factor, A, B, offset, tau_c_ei, v0, v_max, weights_before_redistribution
    real(dp) :: standard_deviation, dummy
    real(dp), dimension(:,:,:), allocatable :: rand_matrix
    character(len=100) :: filename, ns_str
    real(dp), dimension(:), allocatable :: data_for_diffusion_coefficient, data_for_convection_coefficient
    logical :: ignore_condition

    s%n_particles = 40000
    call random_number(dummy)

    if (.not.allocated(rand_matrix)) allocate(rand_matrix(5,s%n_particles,1))
    if (.not.allocated(dc%s_vertices)) allocate(dc%s_vertices(grid_size(1)+1))
    if (.not.allocated(dc%A)) allocate(dc%A(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    if (.not.allocated(dc%B)) allocate(dc%B(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    if (.not.allocated(dc%grad_A)) allocate(dc%grad_A(grid_size(1)))
    if (.not.allocated(dc%grad_B)) allocate(dc%grad_B(grid_size(1)))
    dc%s_vertices = verts_sthetaphi(1, [(grid_size(3)*(ns-1)+1,ns=1,grid_size(1)+1)])
    s%k = 1000
    s%j = 100
    if (.not.allocated(s%delta_s)) allocate(s%delta_s(s%k))
    if (.not.allocated(s%delta_s_squared)) allocate(s%delta_s_squared(s%k))
    if (.not.allocated(s%time)) allocate(s%time(s%k))
    if (.not.allocated(s%f_v)) allocate(s%f_v(s%k,s%j))
    if (.not.allocated(s%check)) allocate(s%check(s%k))
    if (.not.allocated(data_for_diffusion_coefficient)) allocate(data_for_diffusion_coefficient(s%n_particles))
    if (.not.allocated(data_for_convection_coefficient)) allocate(data_for_convection_coefficient(s%n_particles))
    if (.not.allocated(s%boole_large_distance)) allocate(s%boole_large_distance(s%n_particles))
    s%boole_large_distance = .false.

    call allocate_start_type(s%n_particles)
    call set_particle_type_specifications
    start%v0(2) = start%v0(2)*sqrt(s%temperature/in%energy_eV)
    call initialize_exit_data(s%n_particles)

    tau_c_ei = 1.7_dp*1.0d-4 !rough estimate from nrl formula booklet
    start%t(2) = 2.0d0*tau_c_ei !check afterwards if this was too little time

    s%time = [(start%t(2)/s%k*i,i = 1,s%k)]

    ns_start = 15

    do ns = ns_start,grid_size(1)

        print*, 'ns = ', ns, '/', grid_size(1)


        !Initiate electrons at the different flux surfaces leaving out the boundaries
        s%s0 = dc%s_vertices(ns)
        s%s0 = 0.5d0

        call RANDOM_NUMBER(rand_matrix)
        call set_starting_positions(rand_matrix,(/2/), s%s0)
        call set_rest_of_individual_particle_specifications(rand_matrix,boole_diffusion_coefficient_in = .true., &
                                                            species_in=(/2/),n_particles_in=s%n_particles)
        call set_weight


        call prepare_next_round_of_parallelised_particle_pushing(2)

        if (ns.eq.ns_start) then 
            call calc_collision_coefficients_for_all_tetrahedra(2)

            v0 = start%v0(2)
            v_max = v0*sqrt(start%epsilon_max)
            weights_before_redistribution = g%total_volume*in%density
            call prepare_russian_roulette(v0,v_max,weights_before_redistribution,10)
        endif

        ! open(23,file = 'energy_before_pushing.dat')
        !     do i = 1,s%n_particles
        !         write(23,*) start%energy(i,2)
        !     enddo
        ! close(23)
        open(71)
        call parallelised_particle_pushing(species = 2,j = 1,boole_diffusion_coefficient = .true.,n_particles_in=s%n_particles)
        close(71)



        ! open(23,file = 'energy_after_pushing.dat')
        !     do i = 1,s%n_particles
        !         write(23,*) start%energy(i,2)
        !     enddo
        ! close(23)

        !Starting index (Throw away first 20 percent of values)
        i = ceiling(dble(s%k)*0.2_dp)
        

        call llsq ( int(s%k - i + 1, kind=8), s%time(i:s%k), s%delta_s(i:s%k), A, offset)

        data_for_diffusion_coefficient = s%delta_s_squared-2*A*s%time*s%delta_s+A**2*s%time**2
        call llsq ( int(s%k - i + 1, kind=8), s%time(i:s%k), data_for_diffusion_coefficient, B, offset)
        B = B/2

        ! n_ignored = 0
        ! data_for_convection_coefficient = exit_data%x(1,:,2)-s%s0
        ! do i = 1,s%n_particles
        !     ignore_condition = abs(data_for_convection_coefficient(i)).gt.0.99_dp*min(s%s0-sfc_s_min,1.0_dp-s%s0)
        !     if (ignore_condition) then
        !         data_for_convection_coefficient(i) = 0.0_dp
        !         n_ignored = n_ignored + 1
        !     endif
        ! enddo

        ! A = sum(exit_data%x(1,:,2)-s%s0)/(s%n_particles*start%t(2))
        ! data_for_diffusion_coefficient = abs(exit_data%x(1,:,2)-(s%s0+A*start%t(2)))
        ! call quicksort(data_for_diffusion_coefficient, 1, s%n_particles)
        ! standard_deviation = data_for_diffusion_coefficient(int(s%n_particles*0.682689_dp))
        ! B = standard_deviation**2/(2*start%t(2))


        dc%A(ns) = A
        dc%B(ns) = B

        print*, 'A = ', A 
        print*, 'B = ', B
    write(ns_str, '(I0)') ns
    filename = 'exit_s_values' // trim(ns_str) // '.dat'
    open(newunit=file_id,file = filename)
        do i = 1,s%n_particles
            write(file_id,*) exit_data%x(1,i,2)
        enddo
    close(file_id)


    ! filename = 'data' // trim(ns_str) // '.dat'
    ! open(23,file = filename)
    !     do i = 1,s%k
    !         write(23,*) s%time(i), s%delta_s(i), s%delta_s_squared(i)
    !     enddo
    ! close(23)

    ! open(23,file = 'f_v.dat')
    !     do i = 1,s%j
    !         write(23,*) s%f_v(:,i)
    !     enddo
    ! close(23)

    !print*, start%weight(:,2)

    enddo

    !Extrapolate values to the boundaries
    extrapolation_factor = (dc%s_vertices(2)-dc%s_vertices(1))/(dc%s_vertices(3)-dc%s_vertices(2))
    dc%A(1) =     dc%A(2) + (dc%A(2)-dc%A(3))*extrapolation_factor
    dc%B(1) = max(dc%B(2) + (dc%B(2)-dc%B(3))*extrapolation_factor, 0.0_dp)
    extrapolation_factor = (dc%s_vertices(grid_size(1)+1)-dc%s_vertices(grid_size(1)))/ &
                           (dc%s_vertices(grid_size(1))-dc%s_vertices(grid_size(1)-1))
    dc%A(grid_size(1)+1) =     dc%A(grid_size(1)) + (dc%A(grid_size(1))-dc%A(grid_size(1)-1))*extrapolation_factor
    dc%B(grid_size(1)+1) = max(dc%B(grid_size(1)) + (dc%B(grid_size(1))-dc%B(grid_size(1)-1))*extrapolation_factor, 0.0_dp)

    !Calculate gradients
    do ns = 1,grid_size(1)
        dc%grad_A(ns) = (dc%A(ns+1)-dc%A(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
        dc%grad_B(ns) = (dc%B(ns+1)-dc%B(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
    enddo

    open(23,file = 'A_and_B.dat')
    do i = 1,grid_size(1)+1
        write(23,*) dc%A(i), dc%B(i)
    enddo
    close(23)

end subroutine calc_electron_diffusion_coefficients

subroutine calc_electron_density_via_random_walk(iteration_step) !call this after every ion pushing sequence

    use gorilla_applets_types_mod, only: in, time_t, dc, ep, g, output, start, exit_data
    use tetra_grid_settings_mod, only: grid_size, sfc_s_min
    use binsrc_mod, only: binsrc
    use tetra_grid_mod, only: ntetr, verts_sthetaphi

    real(dp) :: delta_x, delta_t, xi, A, B, cell_size, B_fit, A_fit, p
    integer :: i, ns, k, num_steps_min, count_lost_particles, num_particles, particle_multiplication, iteration_step
    real(dp), dimension(:), allocatable :: electron_density, electron_prism_densities
    real(dp), dimension(:), allocatable :: position, weight, exit_time
    type(time_t) :: t
    logical :: boole_lost, boole_use_fit_function = .true.

    particle_multiplication = max(int(1.0d4/in%num_particles),1)
    num_particles = in%num_particles*particle_multiplication

    allocate(electron_density(grid_size(1)), electron_prism_densities(ntetr/3))
    electron_density = 0.0_dp

    !>This should later be an input
    allocate(position(num_particles))
    allocate(weight(num_particles))
    allocate(exit_time(num_particles))

    do i = 1,particle_multiplication
        position(i:num_particles:particle_multiplication) = start%x(1,:,2)
        weight(i:num_particles:particle_multiplication) = start%weight(:,2)
    enddo

    ! call RANDOM_NUMBER(position)
    ! call generate_distribution_one_minus_x2(position)
    ! allocate(dc%A(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    ! allocate(dc%B(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    ! allocate(dc%grad_A(grid_size(1)))
    ! allocate(dc%grad_B(grid_size(1)))
    ! allocate(dc%s_vertices(grid_size(1)+1))
    ! dc%s_vertices = verts_sthetaphi(1, [(grid_size(3)*(ns-1)+1,ns=1,grid_size(1)+1)])
    ! dc%A = 0.0_dp
    ! dc%B = 0.1_dp
    ! dc%grad_A = 0.0_dp
    ! dc%grad_B = 0.0_dp


    if (iteration_step.eq.1) then 
        call get_diffusion_coefficient_data(boole_use_fit_function)
    else
        call calc_convection_coefficient_from_electric_field(boole_use_fit_function)
    endif

    t%step = in%time_step
    num_steps_min = 1000

count_lost_particles = 0
    do i = 1,num_particles
        t%confined = 0.0_dp
        boole_lost = .false.
        do while ((t%confined.lt.t%step).and.(.not.boole_lost))

            !binsrc finds k such that dc%s_vertices(k-1) < position(i) < dc%s_vertices(k)
            call binsrc(dc%s_vertices,1,grid_size(1)+1,position(i),k) 
            k = k-1 
            cell_size = dc%s_vertices(k+1) - dc%s_vertices(k)
            A = dc%A(k) + dc%grad_A(k)*(position(i)-dc%s_vertices(k))
            B = dc%B(k) + dc%grad_B(k)*(position(i)-dc%s_vertices(k))

            if (boole_use_fit_function) then
                B_fit = sum(dc%polynomial_coefficients_for_B*(/1.0_dp,position(i),position(i)**2,0.0_dp/))
                A_fit = sum(dc%polynomial_coefficients_for_A*(/1.0_dp,position(i),0.0_dp/))
                if (position(i).lt.dc%s_vertices(4)) then
                    p = (position(i) - sfc_s_min)/(dc%s_vertices(4) - sfc_s_min)
                    B = p*B_fit + (1-p)*B
                    A = p*A_fit + (1-p)*B/position(i) + A
                else
                    B = B_fit
                    A = A + A_fit
                endif
            endif


            delta_t = min(t%step/num_steps_min, cell_size**2/((abs(A)+abs(B))*4), t%step - t%confined) !control maximum possible jump

            electron_density(k) = electron_density(k) + delta_t*weight(i)
            t%confined = t%confined + delta_t

            call random_number(xi)
            xi = sqrt(12.0_dp)*(xi-0.5_dp)
            delta_x = sqrt(2*delta_t*B)*xi + A*delta_t
            !print*, 'diffusive part over convective part = ', A,B,abs(sqrt(2*delta_t*B)/(A*delta_t)), &
            !sqrt(2*delta_t*B)*xi/(A*delta_t)

            position(i) = position(i) + delta_x

            if (position(i).lt.sfc_s_min) then 
                position(i) = 2.0_dp*sfc_s_min - position(i)
            endif
            if (position(i).gt.1.0_dp)    then
                boole_lost = .true.
                count_lost_particles = count_lost_particles+1
            endif
        enddo
        exit_time(i) = t%confined
    enddo

    if (count_lost_particles.lt.num_particles) print*, 'Warning: the tracing time (', t%step,'s) was so short that only ',&
                                                        count_lost_particles,'out of', num_particles, &
                                                        'electrons left the computation domain'

    open(10)
    do ns = 1, grid_size(1)
        write(10,*) sfc_s_min + (1.0_dp - sfc_s_min) * (ns - 0.5_dp) / (grid_size(1)+1), electron_density(ns)
    enddo
    close(10)

    do i = 1,grid_size(1)
        electron_density(i) = electron_density(i)/(ep%s_shell_volumes(i)*t%step*num_particles)
        call fill_vector_parts_with_value(electron_prism_densities, g%prisms_per_flux_tube(i,:), electron_density(i))
    enddo


    !This is a bit cumbersome having to use electron_prism_densities and a loop over all prisms instead of simply putting 
    !output%prism_moments into fill_vector_parts_with_value, but that does not work because the latter does not accept 
    !a complex argument. Think about a solution to this problem
    do i = 1,ntetr/3
        output%prism_moments(1,i,2) = complex(electron_prism_densities(i), 0.0_dp)
    enddo

    print*, 'Total tracing time of all electrons divided by number of particles is: ', &
    sum(exit_time)/num_particles, 's'
    !print*, 'exit times are : ', exit_data%t_confined(:,2)

end subroutine calc_electron_density_via_random_walk

subroutine get_diffusion_coefficient_data(boole_use_fit_function)

    use gorilla_applets_types_mod, only: dc
    use tetra_grid_settings_mod, only: grid_size
    use tetra_grid_mod, only: verts_sthetaphi
    use utils_polyfit_mod, only: quadratic_fit

    character(len=100) :: filename
    integer :: id, ns, lower_bound, upper_bound
    real(dp), dimension(:), allocatable :: x,y
    logical :: success, boole_use_fit_function
    real(dp) :: a0, a1, a2


    if (.not.allocated(dc%A)) allocate(dc%A(grid_size(1)+1))
    if (.not.allocated(dc%A_from_first_run)) allocate(dc%A_from_first_run(grid_size(1)+1))
    if (.not.allocated(dc%B)) allocate(dc%B(grid_size(1)+1))
    if (.not.allocated(dc%grad_A)) allocate(dc%grad_A(grid_size(1)))
    if (.not.allocated(dc%grad_B)) allocate(dc%grad_B(grid_size(1)))
    if (.not.allocated(dc%s_vertices)) allocate(dc%s_vertices(grid_size(1)+1))

    filename = 'A_and_B.dat'
    open(newunit=id,file = filename,action='read')
    do ns = 1,grid_size(1)+1
        read(id,*) dc%A(ns), dc%B(ns)
    enddo

    dc%s_vertices = verts_sthetaphi(1, [(grid_size(3)*(ns-1)+1,ns=1,grid_size(1)+1)])

    dc%A_from_first_run = dc%A
    if (boole_use_fit_function) dc%A = 0.0_dp

    do ns = 1,grid_size(1)
        dc%grad_A(ns) = (dc%A(ns+1)-dc%A(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
        dc%grad_B(ns) = (dc%B(ns+1)-dc%B(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
    enddo



    lower_bound = 2
    upper_bound = grid_size(1) - 4

    if (.not.allocated(x)) allocate(x(upper_bound-lower_bound+1))
    if (.not.allocated(y)) allocate(y(upper_bound-lower_bound+1))

    x = dc%s_vertices(lower_bound:upper_bound)
    y = dc%B(lower_bound:upper_bound)

    call quadratic_fit(size(x), x, y, a0, a1, a2, success) ! Least-squares fit of y ≈ a2*x^2 + a1*x + a0

    dc%polynomial_coefficients_for_B = (/a0,a1,a2,0.0_dp/)
    dc%polynomial_coefficients_for_A = (/a1,2*a2,0.0_dp/)

end subroutine get_diffusion_coefficient_data

subroutine calc_convection_coefficient_from_electric_field(boole_use_fit_function)

    use tetra_physics_mod, only: tetra_physics
    use constants, only: echarge, ev2erg
    use gorilla_applets_types_mod, only: dc, in, s
    use tetra_grid_settings_mod, only: grid_size

    logical :: boole_use_fit_function
    integer :: i, ns
    real(dp) :: electric_field

    do i = 1,grid_size(1)+1  
        if (i.eq.1) then
            electric_field = -tetra_physics(1)%gPhi(1)
        elseif (i.eq.grid_size(1)+1) then
            electric_field = -tetra_physics((grid_size(1)-1)*6*grid_size(2)+1)%gPhi(1)
        else
            electric_field = -0.5_dp*(tetra_physics((i-2)*6*grid_size(2)+1)%gPhi(1) + &
                                      tetra_physics((i-1)*6*grid_size(2)+1)%gPhi(1))
        endif
        dc%A(i) = -dc%B(i)*echarge*electric_field/(s%temperature*ev2erg)
        if (.not.boole_use_fit_function) dc%A(i) = dc%A(i) + dc%A_from_first_run(i)
    enddo

    do ns = 1,grid_size(1)
        dc%grad_A(ns) = (dc%A(ns+1)-dc%A(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
        dc%grad_B(ns) = (dc%B(ns+1)-dc%B(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
    enddo

end subroutine calc_convection_coefficient_from_electric_field

subroutine generate_distribution_x4_exp_neg_x(b, output_array)

     use gorilla_applets_types_mod, only: s, in

    real(dp), dimension(:,:), intent(out) :: output_array
    real(dp), intent(in) :: b
    
    real(dp), dimension(:), allocatable :: uniform_random
    real(dp) :: x, y, max_value
    integer :: j, accept_count, i, num_cols, row, num_rows

    num_cols = size(output_array, 1)
    num_rows = size(output_array, 2)
    allocate(uniform_random(num_cols))

    max_value = 4.0_dp**4.0_dp * (s%temperature/in%energy_eV)**4.0_dp*exp(-4.0_dp)
    
    do row = 1, num_rows
        accept_count = 0
        
        do while (accept_count < num_cols)
            call random_number(uniform_random)
            
            do j = 1, num_cols
                if (accept_count >= num_cols) exit
                
                x = b * uniform_random(j)
                y = x**4.0_dp * exp(-x/(s%temperature/in%energy_eV))
                
                call random_number(uniform_random(j))
                if (uniform_random(j) * max_value <= y) then
                    accept_count = accept_count + 1
                    output_array(accept_count, row) = x
                endif
            enddo
        enddo
    enddo

end subroutine generate_distribution_x4_exp_neg_x

subroutine generate_distribution_sqrt_x_exp_neg_x(xmax, output_array)

    real(dp), dimension(:,:), intent(out) :: output_array
    real(dp), intent(in) :: xmax
    
    real(dp), dimension(:), allocatable :: uniform_random
    real(dp) :: x, y, max_value
    integer :: j, accept_count, i, num_cols, row, num_rows

    num_cols = size(output_array, 1)
    num_rows = size(output_array, 2)
    allocate(uniform_random(num_cols))

    max_value = sqrt(0.5_dp) * exp(-0.5_dp)
    
    do row = 1, num_rows
        accept_count = 0
        
        do while (accept_count < num_cols)
            call random_number(uniform_random)
            
            do j = 1, num_cols
                if (accept_count >= num_cols) exit
                
                x = xmax * uniform_random(j)
                y = sqrt(x) * exp(-x)
                
                call random_number(uniform_random(j))
                if (uniform_random(j) * max_value <= y) then
                    accept_count = accept_count + 1
                    output_array(accept_count, row) = x
                endif
            enddo
        enddo
    enddo

end subroutine generate_distribution_sqrt_x_exp_neg_x

subroutine generate_marker_distribution(xmax, output_array)

    use binsrc_mod, only: binsrc

!for a velocity distribution according to (1+v^7)*v^2*exp(-v^2),
!we have an energy distribution according to (1+E^(7/2))*sqrt(E)*exp(-E)

    real(dp), dimension(:,:), intent(out) :: output_array
    real(dp), intent(in) :: xmax
    integer, parameter :: nsize=10000
    logical, save :: firstentry = .true.
    real(dp), dimension(0:nsize), save :: pdf
    integer :: i,j, k, num_cols, num_rows
    real(dp) :: x,xi
    real(dp), save :: hx
    
    if(firstentry) then
        hx=xmax/dble(nsize)
        firstentry = .false.
        pdf(0)=0.d0
        do i=1,nsize
            x=(dble(i)-0.5d0)*hx
            pdf(i)=pdf(i-1)+(1.d0+x**7)*x**2*exp(-x**2)
        enddo
        pdf=pdf/pdf(nsize)
    endif

    num_cols = size(output_array, 1)
    num_rows = size(output_array, 2)
    
    do i = 1,num_cols
        do j = 1, num_rows
            call random_number(xi)
            call binsrc(pdf,0,nsize,xi,k)
            output_array(i,j) = (dble(k)-0.5d0)*hx
        enddo
    enddo

end subroutine generate_marker_distribution

subroutine generate_distribution_one_minus_x2(output_array)

    use binsrc_mod, only: binsrc

    real(dp), dimension(:,:), intent(out) :: output_array

    integer, parameter :: nsize = 10000
    logical, save :: firstentry = .true.
    real(dp), dimension(0:nsize), save :: cdf
    real(dp), save :: hx

    integer :: i, j, k, num_cols, num_rows
    real(dp) :: xi, x, fx

    if (firstentry) then
        hx = 1.0_dp / dble(nsize)
        firstentry = .false.
        cdf(0) = 0.0_dp
        do i = 1, nsize
            x  = (dble(i) - 0.5_dp) * hx
            fx = 1.0_dp - x*x
            cdf(i) = cdf(i-1) + fx
        enddo
        cdf = cdf / cdf(nsize)
    endif

    num_cols = size(output_array, 1)
    num_rows = size(output_array, 2)

    do i = 1, num_cols
        do j = 1, num_rows
            call random_number(xi)
            call binsrc(cdf, 0, nsize, xi, k)
            output_array(i, j) = (dble(k) - 0.5_dp) * hx
        enddo
    enddo

end subroutine generate_distribution_one_minus_x2

subroutine sort_array(array)
  real(dp), dimension(:) :: array
  integer, parameter :: n = 6
  real :: a(n) = [3.5, 1.2, -4.0, 7.8, 0.0, 2.1]
  integer :: i, j
  real(dp) :: temp

  ! simple bubble sort (ascending)
  do i = 1, size(array)
     do j = i+1, size(array)
        if (array(j) < array(i)) then
           temp = array(i)
           array(i) = array(j)
           array(j) = temp
        end if
     end do
  end do

end subroutine sort_array

recursive subroutine quicksort(arr, left, right)
real(dp), intent(inout) :: arr(:)
integer, intent(in) :: left, right
integer :: i, j
real(dp) :: pivot, temp

i = left
j = right
pivot = arr((left+right)/2)

do
    do while (arr(i) < pivot)
        i = i + 1
    end do
    do while (arr(j) > pivot)
        j = j - 1
    end do
    if (i <= j) then
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
        i = i + 1
        j = j - 1
    end if
    if (i > j) exit
end do

if (left < j) call quicksort(arr, left, j)
if (i < right) call quicksort(arr, i, right)
end subroutine quicksort

end module utils_self_consistent_ef_mod
