module utils_scef_particle_pushing_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine parallelised_particle_pushing(species,j,boole_diffusion_coefficient,n_particles_in)

    use gorilla_applets_types_mod, only: counter, c, in, time_t, moment_specs, counter_t, particle_status_t, start, exit_data, s, &
    g, maximum_s, weights
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
    !$OMP& SHARED(counter, kpart,species, in, c, iantithetic, start, j, s, boole_diffusion_coefficient,n_particles,rr,weights) &
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
                        call play_russian_roulette(weights%w(n,species),v,v_save,particle_state_for_rr,local_rr)
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
    use gorilla_applets_types_mod, only: counter_t, particle_status_t, g, start, s, in, time_t, maximum_s, weights
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
        if ((.not.in%boole_static_ne).and.(species.eq.2)) then
            start%x(:,n,2) = start%x(:,n,1)
            weights%w(n,2) = weights%w(n,1)
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
                s%delta_s(k) = s%delta_s(k) + (x(1) - s%s0)*weights%w(n,species)
                s%delta_s_squared(k) = s%delta_s_squared(k) + (x(1) - s%s0)**2*weights%w(n,species)
                if (s%boole_large_distance(n)) then !delete contribution to delta_s, double contribution to delta_s_squared
                    s%delta_s(k) = s%delta_s(k) - (x(1) - s%s0)*weights%w(n,species)
                    s%delta_s_squared(k) = s%delta_s_squared(k) + (x(1) - s%s0)**2*weights%w(n,species)
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
            call play_russian_roulette(weights%w(n,species),v,v_save,particle_state_for_rr,local_rr)
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

    use gorilla_applets_types_mod, only : particle_status_t, time_t, weights
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
        weights%w(n,species) = local_rr%weight(id)

        particle_status%lost = .false.
        particle_status%initialized = .true.
        particle_status%exit = .false.
    endif

end subroutine initiate_next_split_particle

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

subroutine treat_particles_that_are_lost_but_should_not_be(z_save_at_x_save,ind_tetr_save,z_save,x_save,x,vpar,vperp, &
                                                            perpinv,ind_tetr,vpar_save,vperp_save,vpar_init,vperp_init)

    use gorilla_applets_types_mod, only: filenames
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
    open(newunit = problem_unit, file = filenames%pushing_problems, position = 'append')
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

    use gorilla_applets_types_mod, only: in, start, s, weights
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
    weights%w(n,species) = weights%w(n,species)*abs((tetra_physics(ind_tetr)%sqg1 + sum(tetra_physics(ind_tetr)%gsqg*z_save)))
    !print*, 'weight before = ', weights%w(n,species), n

    if (in%boole_linear_density_simulation) then
        weights%w(n,species) = weights%w(n,species)*(1.0_dp-0.9_dp*x(1))
    endif

    if (boole_diffusion_coefficient) then
        weights%w(n,species) = 1.0_dp/s%n_particles
    elseif (.not. in%boole_monoenergetic) then

        temperature = in%energy_eV
        if (boole_diffusion_coefficient) temperature = s%temperature
        weights%w(n,species) = weights%w(n,species)/  &
        (sqrt(start%energy(n,species)/temperature)*exp(-start%energy(n,species)/temperature)/(temperature*ev2erg*sqrt(pi)/2))
        !the last term is the integral of the function from zero to inf over energy in correct units
        !weights%w(n,species) = weights%w(n,species)/((1.0_dp+(start%energy(n,species)/s%temperature)**(3.5_dp))*  &
        !sqrt(start%energy(n,species)/s%temperature)*exp(-start%energy(n,species)/s%temperature)/(s%temperature*ev2erg*14035.0_dp))
        !the last term is the integral of the function from zero to inf over energy in correct units

        phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z_save)
        phi_elec_func = 0.0_dp !when working with fixed sources, electric potentials are not useful since they change the weights
        !and thus the magnitude o the sources

        weights%w(n,species) = weights%w(n,species)*2/sqrt(pi)*sqrt(start%energy(n,species)*ev2erg)
        ! if (.not.boole_diffusion_coefficient) then
        !     weights%w(n,species) = weights%w(n,species)*start%epsilon_max*in%energy_eV*ev2erg
        ! endif
        if (.not. in%boole_linear_temperature_simulation) then
            temperature = in%energy_eV*ev2erg
            if (boole_diffusion_coefficient) temperature = s%temperature*ev2erg
            weights%w(n,species) =weights%w(n,species)/temperature**1.5_dp* &
                        & exp(-(start%energy(n,species)*ev2erg+start%particle_charge(species)*phi_elec_func)/temperature)
        else
            temperature = in%energy_eV*ev2erg*(1.0_dp-0.9*x(1))
            if (boole_diffusion_coefficient) temperature = s%temperature*ev2erg*(1.0_dp-0.9*x(1))
            weights%w(n,species) = weights%w(n,species)/temperature**1.5_dp* &
                        & exp(-(start%energy(n,species)*ev2erg+start%particle_charge(species)*phi_elec_func)/temperature)
        endif
    endif

    start%jperp(n,species) = start%particle_mass(species)*vperp**2*start%cm_over_e(species)/(2*bmod_func(z_save,ind_tetr))*(-1)


    !-1 because of negative gyrophase

!print*, 'weight after = ', weights%w(n,species)*(1.0_dp+(start%energy(n,species)/s%temperature)**(3.5_dp))/14035.0_dp/2*sqrt(pi)
!print*, 'weight after for real = ', weights%w(n,species)



end subroutine calc_particle_weights_and_jperp

end module utils_scef_particle_pushing_mod
