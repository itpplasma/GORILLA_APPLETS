module utils_parallelised_particle_pushing_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64

    implicit none

    !Deferred background updates of the energy and momentum conserving collisions: while
    !boole_deferred_background_updates is set, collisions read a frozen background and every thread
    !adds the parallel momentum and energy [erg] its markers gain into its own slot
    !background_exchange(1:2, species, tetrahedron, thread). apply_deferred_background_updates sums the
    !slots in thread order (reproducible for a fixed thread count) and updates each reservoir once.
    logical :: boole_deferred_background_updates = .false.
    real(dp), dimension(:,:,:,:), allocatable :: background_exchange

contains

subroutine print_progress(num_particles,kpart,n)

    integer, parameter :: n_critical = 100
    integer :: num_particles, kpart, n
    logical :: print_progress_for_every_particle = .false.

    if ((.not.print_progress_for_every_particle).and.(num_particles.gt.n_critical)) then
        if (modulo(kpart,int(num_particles/n_critical)).eq.0) then
            print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
        endif
    else
        print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
    endif

end subroutine print_progress

subroutine handle_lost_particles(local_counter, boole_particle_lost)

    use gorilla_applets_types_mod, only: counter_t

    type(counter_t) :: local_counter
    logical :: boole_particle_lost

    !write another if clause (if hole size = minimal .and. particle lost inside .and. boole_cut_out_hole = .true. 
    !(use extra variable m in orbit routine (0 normally, 1 when lost outside, -1 when lost inside))),
    !if m = 1 do as now, if m = -1 select arbitrary newp position and update x, vpar and vperp)

    !$omp critical
    local_counter%lost_particles = 1
    boole_particle_lost = .true.
    !$omp end critical

end subroutine handle_lost_particles

subroutine initialise_loop_variables(l, n, local_counter,particle_status,t,local_tetr_moments,x,vpar,vperp,species_in)

use gorilla_applets_types_mod, only: in, counter_t, start, time_t, particle_status_t
use constants, only: ev2erg
use utils_orbit_timestep_mod, only: is_touched
use gorilla_applets_settings_mod, only: i_option

integer, intent(in) :: l, n
integer, intent(in), optional :: species_in
type(counter_t) :: local_counter
type(particle_status_t) :: particle_status
type(time_t) :: t
real(dp) :: t_step, t_confined, pitchpar, v, vpar, vperp, epsilon
real(dp), dimension(3) :: x
complex(dp), dimension(:,:) :: local_tetr_moments
integer :: species
species = 1
if (present(species_in)) species = species_in

call set_counter_zero(local_counter)
particle_status%lost = .false.
particle_status%initialized = .false.
particle_status%exit = .false.
t%step = in%time_step
if (present(species_in)) t%step = start%t(species)
t%confined = 0.0_dp
!with touched-tetrahedron tracking, add_marker_moments_to_thread_sums resets the moments instead
if ((l.eq.1).and.(.not.allocated(is_touched))) local_tetr_moments = 0.0_dp
pitchpar = start%pitch(n,species)
x = start%x(:,n,species)
v = sqrt(start%energy(n,species)*ev2erg*2/start%particle_mass(species))
vpar = pitchpar * v
vperp = sqrt(v**2-vpar**2)

epsilon = 1.0d-1 !ensures that in one step, a maximum displacement of epsilon cm occurs
if (i_option.eq.13) then
    if (in%anomalous_diffusion_coefficient > 0.0d0) then
        t%step_anomalous_transport = min(epsilon**2/(2.0d0*in%anomalous_diffusion_coefficient),t%step/300.0d0)
    else
        t%step_anomalous_transport = t%step  ! No anomalous transport, use full time step
    endif
endif

end subroutine initialise_loop_variables

subroutine add_marker_moments_to_thread_sums(local_tetr_moments, thread_tetr_moments, thread_prism_moments_squared)

    !Adds the moments of one marker (or antithetic pair) to the sums of the calling thread, without
    !synchronisation. The squared moments need the square of the marker's total per prism (its
    !contributions to a prism are correlated), so they are formed here, marker by marker.
    !With touched-tetrahedron tracking, only the visited tetrahedra are transferred and then reset.

    use gorilla_applets_types_mod, only: moment_specs
    use utils_orbit_timestep_mod, only: touched_tetr, is_touched, n_touched

    complex(dp), dimension(:,:), intent(inout) :: local_tetr_moments
    complex(dp), dimension(:,:), intent(inout) :: thread_tetr_moments, thread_prism_moments_squared
    integer :: i, k, ind_tetr, first_in_prism

    if (.not.allocated(is_touched)) then
        thread_tetr_moments = thread_tetr_moments + local_tetr_moments
        if (moment_specs%boole_squared_moments) then
            do k = 1, size(thread_prism_moments_squared, 2)
                thread_prism_moments_squared(:,k) = thread_prism_moments_squared(:,k) + &
                    (local_tetr_moments(:,3*k-2) + local_tetr_moments(:,3*k-1) + local_tetr_moments(:,3*k))**2
            enddo
        endif
        return
    endif

    do i = 1, n_touched
        ind_tetr = touched_tetr(i)
        thread_tetr_moments(:,ind_tetr) = thread_tetr_moments(:,ind_tetr) + local_tetr_moments(:,ind_tetr)
        if (moment_specs%boole_squared_moments) then
            !square each touched prism once, at its first touched tetrahedron
            k = (ind_tetr + 2)/3
            first_in_prism = 3*k - 2
            do while (.not.is_touched(first_in_prism))
                first_in_prism = first_in_prism + 1
            enddo
            if (first_in_prism.eq.ind_tetr) thread_prism_moments_squared(:,k) = thread_prism_moments_squared(:,k) + &
                (local_tetr_moments(:,3*k-2) + local_tetr_moments(:,3*k-1) + local_tetr_moments(:,3*k))**2
        endif
    enddo
    do i = 1, n_touched
        local_tetr_moments(:,touched_tetr(i)) = (0.0_dp, 0.0_dp)
        is_touched(touched_tetr(i)) = .false.
    enddo
    n_touched = 0

end subroutine add_marker_moments_to_thread_sums

subroutine add_thread_moments_to_output(thread_tetr_moments, thread_prism_moments_squared, species)

    !Adds the moment sums of one thread to the output; call once per thread inside a critical section.

    use gorilla_applets_types_mod, only: moment_specs, output

    complex(dp), dimension(:,:), intent(in) :: thread_tetr_moments, thread_prism_moments_squared
    integer, intent(in) :: species
    integer :: k

    output%tetr_moments(:,:,species) = output%tetr_moments(:,:,species) + thread_tetr_moments
    do k = 1, size(output%prism_moments, 2)
        output%prism_moments(:,k,species) = output%prism_moments(:,k,species) + &
            thread_tetr_moments(:,3*k-2) + thread_tetr_moments(:,3*k-1) + thread_tetr_moments(:,3*k)
    enddo
    if (moment_specs%boole_squared_moments) then
        output%prism_moments_squared(:,:,species) = output%prism_moments_squared(:,:,species) + thread_prism_moments_squared
    endif

end subroutine add_thread_moments_to_output

subroutine set_counter_zero(counter)

    use gorilla_applets_types_mod, only: counter_t
    
    type(counter_t), intent(inout) :: counter
    
    counter%lost_particles = 0
    counter%lost_inside = 0
    counter%tetr_pushings = 0
    counter%phi_0_mappings = 0
    counter%integration_steps = 0
    counter%pusher_failures = 0

end subroutine set_counter_zero

subroutine add_local_counter_to_counter(local_counter)

    use gorilla_applets_types_mod, only: counter_t, counter
    
    type(counter_t), intent(in) :: local_counter
    
    counter%lost_particles = counter%lost_particles + local_counter%lost_particles
    counter%lost_inside = counter%lost_inside + local_counter%lost_inside
    counter%tetr_pushings = counter%tetr_pushings + local_counter%tetr_pushings
    counter%phi_0_mappings = counter%phi_0_mappings + local_counter%phi_0_mappings
    counter%pusher_failures = counter%pusher_failures + local_counter%pusher_failures
    
end subroutine add_local_counter_to_counter

subroutine carry_out_collisions(i, n, t, x, vpar, vperp, ind_tetr, iface, species_in, iswmode_in)

    use gorilla_applets_types_mod, only: in, time_t
    use find_tetra_mod, only: find_tetra

    integer, intent(in) :: i, n
    integer, intent(in), optional :: species_in
    integer, intent(in), optional :: iswmode_in
    integer :: species
    integer :: iswmode
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(inout) :: vpar, vperp
    type(time_t) :: t
    integer :: ind_tetr, iface

    species = 1

    if (present(species_in)) species = species_in
    !iswmode options:
    !1 - full operator (pitch-angle and energy scattering and drag)
    !2 - energy scattering and drag only
    !3 - drag only
    !4 - pitch-angle scattering only
    iswmode = 1
    if (present(iswmode_in)) iswmode = iswmode_in

    if (i.eq.1) call find_tetra(x,vpar,vperp,ind_tetr,iface)
    if (.not.(ind_tetr.eq.-1)) then
        if (in%boole_preserve_energy_and_momentum_during_collisions) then
            call collisions_with_background_updates(i, n, t, x, vpar, vperp, ind_tetr, species, iswmode)
        else
            call collisions_without_background_updates(i, n, t, x, vpar, vperp, ind_tetr, species, iswmode)
        endif
    endif

end subroutine carry_out_collisions

subroutine collisions_with_background_updates(i, n, t, x, vpar, vperp, ind_tetr, species, iswmode)

    use gorilla_applets_types_mod, only: in, c, time_t, start, weights
    use collis_ions, only: stost
    use collis_ions, only: collis_init
    use tetra_physics_mod, only: particle_mass,particle_charge
    use constants, only: echarge,amp
    use gorilla_applets_settings_mod, only: i_option
    use omp_lib, only: omp_get_thread_num

    integer, intent(in) :: i, n, species, iswmode
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(inout) :: vpar, vperp
    type(time_t) :: t
    integer :: ind_tetr

    real(dp), dimension(5) :: zet
    real(dp), dimension(3) :: randnum
    real(dp), dimension(1) :: m, z, dens, temp, efcolf,velrat,enrat
    real(dp) :: vpar_background
    real(dp) :: m0, z0, vpar_save, vperp_save, delta_epsilon, delta_vpar
    integer :: err, j, p, thread
    real(dp) :: marker_to_background_ratio, t_max, w_phys

    do j = 1,c%n-1
        
        m0 = particle_mass
        z0 = particle_charge/echarge
        m = c%mass(j)
        z = c%charge_num(j)
        dens = c%dens_mat(j,ind_tetr)
        temp = c%temp_mat(j,ind_tetr)
        call collis_init(m0,z0, m,z, dens, temp, in%energy_eV, start%v0(species), efcolf, velrat, enrat)

        vpar_save = vpar
        vperp_save = vperp
        vpar_background = c%vpar_mat(j,ind_tetr)
        vpar = vpar - vpar_background

        zet(1:3) = x !spatial position
        zet(4) = sqrt(vpar**2+vperp**2)/start%v0(species) !normalized velocity module 
        zet(5) = vpar/sqrt(vpar**2+vperp**2) !pitch parameter

        t_max = (start%t(species)-t%confined)*start%v0(species)
        if (i_option.eq.13) t_max = min(t_max, t%step_anomalous_transport*start%v0(species))
        
        if (in%boole_precalc_collisions) then
            randnum = c%randcol(n,mod(i-1,c%randcoli)+1,:,species) 
            call stost(efcolf,velrat,enrat,zet,t%step,iswmode,err,t_max,randnum)
        else
            call stost(efcolf,velrat,enrat,zet,t%step,iswmode,err,t_max)
        endif
        
        vpar = zet(5)*zet(4)*start%v0(species)+vpar_background
        vperp = sqrt(1-zet(5)**2)*zet(4)*start%v0(species)
        
        !optionally still change particle_mass, particle_charge and cm_over_e, e.g.:
        !particle_charge = particle_charge + echarge
        !particle_mass = particle_mass + ame - amp
        !cm_over_e = clight*particle_mass/particle_charge

        delta_vpar = vpar - vpar_save
        delta_epsilon = particle_mass/2*(vpar**2 + vperp**2 - vpar_save**2 - vperp_save**2)

        !A marker represents weights%w/n_particles physical particles (the moments are normalised the
        !same way), the reservoir holds dens_mat*tetr_volume background particles
        w_phys = real(weights%w(n,species), dp)/in%n_particles

        if (boole_deferred_background_updates) then
            thread = omp_get_thread_num()
            background_exchange(1,j,ind_tetr,thread) = background_exchange(1,j,ind_tetr,thread) + &
                                                       w_phys*particle_mass*delta_vpar
            background_exchange(2,j,ind_tetr,thread) = background_exchange(2,j,ind_tetr,thread) + w_phys*delta_epsilon
        else
            marker_to_background_ratio = w_phys/(c%dens_mat(j,ind_tetr)*c%tetr_volume(ind_tetr))
            !$omp critical
            call update_background_reservoir(particle_mass, c%mass(j), marker_to_background_ratio, &
                delta_vpar, delta_epsilon, c%vpar_mat(j,ind_tetr), c%temp_mat(j,ind_tetr))
            !$omp end critical
        endif
    enddo

end subroutine collisions_with_background_updates

subroutine start_deferred_background_updates

    !Switches the conserving collisions to deferred background updates (see module header)

    use gorilla_applets_types_mod, only: c
    use tetra_grid_mod, only: ntetr
    use omp_lib, only: omp_get_max_threads

    if (allocated(background_exchange)) deallocate(background_exchange)
    allocate(background_exchange(2, c%n, ntetr, 0:omp_get_max_threads()-1))
    background_exchange = 0.0_dp
    boole_deferred_background_updates = .true.

end subroutine start_deferred_background_updates

subroutine apply_deferred_background_updates(n_rejected)

    !Gives each reservoir cell (c%cell_start/c%cell_members, see build_reservoir_cells) the summed
    !momentum and energy that the markers lost in its tetrahedra since the last call, copies the new
    !state to all its tetrahedra and resets the exchange. n_rejected counts cells whose update would
    !have driven the temperature non-positive and was therefore rejected (conservation is broken there).

    use gorilla_applets_types_mod, only: c
    use tetra_physics_mod, only: particle_mass

    integer, intent(out) :: n_rejected
    real(dp) :: momentum, energy, particles, u, temp, u_old, temp_old
    integer :: j, k, m, i, thread

    n_rejected = 0
    do k = 1, size(c%cell_start) - 1
        do j = 1, c%n-1
            momentum = 0.0_dp
            energy = 0.0_dp
            particles = 0.0_dp
            do m = c%cell_start(k), c%cell_start(k+1) - 1
                i = c%cell_members(m)
                do thread = 0, size(background_exchange, 4) - 1
                    momentum = momentum + background_exchange(1,j,i,thread)
                    energy = energy + background_exchange(2,j,i,thread)
                enddo
                particles = particles + c%dens_mat(j,i)*c%tetr_volume(i)
            enddo
            if ((momentum.eq.0.0_dp).and.(energy.eq.0.0_dp)) cycle
            i = c%cell_members(c%cell_start(k))
            u_old = c%vpar_mat(j,i)
            temp_old = c%temp_mat(j,i)
            u = u_old
            temp = temp_old
            !one marker of unit ratio carrying the totals: r*m_t*delta_vpar = momentum/N_b, r*delta_epsilon = energy/N_b
            call update_background_reservoir(particle_mass, c%mass(j), 1.0_dp/particles, &
                momentum/particle_mass, energy, u, temp)
            if ((u.eq.u_old).and.(temp.eq.temp_old)) n_rejected = n_rejected + 1
            c%vpar_mat(j,c%cell_members(c%cell_start(k):c%cell_start(k+1)-1)) = u
            c%temp_mat(j,c%cell_members(c%cell_start(k):c%cell_start(k+1)-1)) = temp
        enddo
    enddo
    background_exchange = 0.0_dp

end subroutine apply_deferred_background_updates

pure subroutine update_background_reservoir(m_t, m_b, r, delta_vpar, delta_epsilon, u_b, temp_b)
    !Background flow u_b and temperature temp_b [eV] take up the parallel momentum and energy [erg]
    !that a test marker of mass m_t and signed marker-to-background ratio r gains in a collision.

    use constants, only: ev2erg

    real(dp), intent(in) :: m_t, m_b, r, delta_vpar, delta_epsilon
    real(dp), intent(inout) :: u_b, temp_b
    real(dp) :: u_b_old, temp_b_old

    u_b_old = u_b
    temp_b_old = temp_b
    u_b = u_b_old - r*m_t*delta_vpar/m_b
    temp_b = temp_b - (2.0_dp*r*delta_epsilon + m_b*(u_b**2 - u_b_old**2))/(3.0_dp*ev2erg)
    ! Guard against a non-positive background temperature: an exchange that
    ! would drive T_b <= 0 is unphysical, so the reservoir keeps its
    ! previous state instead of propagating an invalid temperature.
    if (temp_b <= 0.0_dp) then
        u_b = u_b_old
        temp_b = temp_b_old
    end if

end subroutine update_background_reservoir

subroutine collisions_without_background_updates(i, n, t, x, vpar, vperp, ind_tetr, species, iswmode)

    use gorilla_applets_types_mod, only: in, c, time_t, start
    use collis_ions, only: stost
    use gorilla_applets_settings_mod, only: i_option
    use profile_data_mod, only: kim_nu_loaded, eval_kim_nu, eval_s_from_psi_pol
    use tetra_physics_mod, only: tetra_physics

    integer, intent(in) :: i, n, species, iswmode
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(inout) :: vpar, vperp
    type(time_t) :: t
    integer :: ind_tetr

    real(dp), dimension(5) :: zet
    real(dp), dimension(3) :: randnum
    real(dp), dimension(c%n) :: efcolf,velrat,enrat,vpar_background
    integer :: err
    real(dp) :: t_max

    real(dp) :: z_cell(3), psi_pol_loc, s_loc, nu_kim_phys, nu_kim_code
    logical  :: use_kim_nu_now

    efcolf = c%efcolf_mat(:,ind_tetr)
    velrat = c%velrat_mat(:,ind_tetr)
    enrat = c%enrat_mat(:,ind_tetr)

    vpar_background = c%vpar_mat(:,ind_tetr)

    vpar = vpar - vpar_background(1) ! TODO(#23): handle multi-species vpar_background
    !since vpar_background actually has num_background_particles entries, consider giving it as an extra
    !optional input variable to stost, before randnum
    zet(1:3) = x !spatial position
    zet(4) = sqrt(vpar**2+vperp**2)/start%v0(species) !normalized velocity module
    zet(5) = vpar/sqrt(vpar**2+vperp**2) !pitch parameter

    t_max = (start%t(species)-t%confined)*start%v0(species)
    if (i_option.eq.13) t_max = min(t_max, t%step_anomalous_transport*start%v0(species))

    ! KIM-nu override: only meaningful for the OU operator (iswmode=5).
    ! Convert the splined nu_e(s) [1/s] to stost code units [1/length]
    ! by dividing by the marker's normalisation speed v0.
    use_kim_nu_now = kim_nu_loaded .and. (iswmode == 5)
    if (use_kim_nu_now) then
        z_cell      = x - tetra_physics(ind_tetr)%x1
        psi_pol_loc = tetra_physics(ind_tetr)%Aphi1 &
                    + sum(tetra_physics(ind_tetr)%gAphi * z_cell)
        s_loc       = eval_s_from_psi_pol(psi_pol_loc)
        nu_kim_phys = eval_kim_nu(s_loc)
        nu_kim_code = nu_kim_phys / start%v0(species)
    end if

    if (in%boole_precalc_collisions) then
        randnum = c%randcol(n,mod(i-1,c%randcoli)+1,:,species)
        if (use_kim_nu_now) then
            call stost(efcolf,velrat,enrat,zet,t%step,iswmode,err,t_max,randnum,nu_override=nu_kim_code)
        else
            call stost(efcolf,velrat,enrat,zet,t%step,iswmode,err,t_max,randnum)
        end if
    else
        if (use_kim_nu_now) then
            call stost(efcolf,velrat,enrat,zet,t%step,iswmode,err,t_max,nu_override=nu_kim_code)
        else
            call stost(efcolf,velrat,enrat,zet,t%step,iswmode,err,t_max)
        end if
    endif

    vpar = zet(5)*zet(4)*start%v0(species)+vpar_background(1)
    vperp = sqrt(1-zet(5)**2)*zet(4)*start%v0(species)

    !optionally still change particle_mass, particle_charge and cm_over_e, e.g.:
    !particle_charge = particle_charge + echarge
    !particle_mass = particle_mass + ame - amp
    !cm_over_e = clight*particle_mass/particle_charge

end subroutine collisions_without_background_updates

subroutine update_exit_data(boole_particle_lost,t_confined,x,vpar,vperp,i,n,phi_0_mappings,species_in,ind_tetr)

    use gorilla_applets_types_mod, only: exit_data, in
    use utils_data_pre_and_post_processing_mod, only: normalised_poloidal_flux
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind

    integer, intent(in)                 :: i, n
    integer, intent(in), optional       :: phi_0_mappings, species_in, ind_tetr
    real(dp), dimension(3), intent(in)  :: x
    real(dp), intent(in)                :: t_confined, vpar, vperp
    logical, intent(in)                 :: boole_particle_lost
    integer                             :: species
    real(dp), dimension(3)              :: z_local
    real(dp)                            :: local_poloidal_flux

    species = 1

    if (present(species_in)) species = species_in

    if (boole_particle_lost.eqv..true.) exit_data%lost(n,species) = 1
    if (boole_particle_lost.eqv..false.) exit_data%lost(n,species) = 0
    exit_data%t_confined(n,species) = t_confined
    exit_data%x(:,n,species) = x
    exit_data%vpar(n,species) = vpar
    exit_data%vperp(n,species) = vperp
    exit_data%integration_step(n,species) = i
    if(present(phi_0_mappings)) exit_data%phi_0_mappings(n,species) = phi_0_mappings

    ! Compute flux surface label (s-coordinate, ranging from sfc_s_min to 1 in flux coordinates)
    if (present(ind_tetr) .and. ind_tetr /= -1) then
        z_local = x - tetra_physics(ind_tetr)%x1
        if (coord_system == 2) then
            ! Flux coordinates: use s-coordinate directly (x(1) is s)
            exit_data%flux_surface(n,species) = tetra_physics(ind_tetr)%x1(1) + z_local(1)
        else if (grid_kind /= 3) then
            ! Cylindrical coordinates with axisymmetric device: use poloidal flux from A_phi
            local_poloidal_flux = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi * z_local)
            exit_data%flux_surface(n,species) = normalised_poloidal_flux(local_poloidal_flux)
        else
            ! grid_kind == 3 (stellarator) with cylindrical coordinates: not supported
            print*, 'Error in update_exit_data: Computing flux surface label from A_phi is only valid for &
                    &axisymmetric devices. For stellarators (grid_kind=3), use flux coordinates (coord_system=2).'
            stop
        endif
    else
        exit_data%flux_surface(n,species) = -1.0_dp  ! Mark as invalid if outside domain
    endif

end subroutine update_exit_data

subroutine update_start_type(x,vpar,vperp,n,species,ind_tetr)

    use gorilla_applets_types_mod, only: start, in
    use supporting_functions_mod, only: bmod_func
    use constants, only: ev2erg
    use tetra_physics_mod, only: tetra_physics

    real(dp), dimension(3), intent(in) :: x
    real(dp), intent(in)               :: vpar, vperp
    integer, intent(in)                :: n, species, ind_tetr
    real(dp), dimension(3)             :: x_local
    real(dp)                           :: v

    start%x(:,n,species) = x
    v = sqrt(vpar**2+vperp**2)
    start%pitch(n,species) = vpar/v
    start%energy(n,species) = 0.5_dp*v**2*start%particle_mass(species)/ev2erg

    if (ind_tetr.eq.-1) then
        start%lost(n,species) = .true.
    else
        x_local = x-tetra_physics(ind_tetr)%x1
        start%jperp(n,species) = start%particle_mass(species)*vperp**2*start%cm_over_e(species)/(2*bmod_func(x_local,ind_tetr))*(-1)
        !-1 because of negative gyrophase
    endif

end subroutine update_start_type

function linspace(start, stop, n) result(x)
    real(dp), intent(in) :: start, stop
    integer, intent(in) :: n
    real(dp) :: x(n)
    real(dp) :: dx
    integer :: i

    dx = (stop - start) / (n - 1)
    x(1) = start
    do i = 2, n
        x(i) = x(i-1) + dx
    end do
end function linspace

subroutine initialise_seed_for_random_numbers_for_each_thread(thread_num, second_factor)
    !Sets an individual, reproducible RNG state for each OpenMP thread. The state
    !is derived from the master seed loaded from seed.inp (via
    !set_seed_for_random_numbers) by combining it with thread_num and an optional
    !second_factor through the SplitMix64 bit mixer, which produces
    !well-decorrelated 64-bit words for neighbouring inputs. Same seed.inp +
    !same thread partition => bit-identical per-thread streams. If no master
    !seed was loaded, we fall back to a hardcoded constant so callers that
    !skipped set_seed_for_random_numbers keep working.

    use utils_data_pre_and_post_processing_mod, only: master_seed_u64, master_seed_loaded, splitmix64_next

    integer, intent(in) :: thread_num
    integer, intent(in), optional :: second_factor
    integer(int64), parameter :: gamma        = -7046029254386353131_int64  ! 0x9E3779B97F4A7C15
    integer(int64), parameter :: second_salt  = -6534898150910275811_int64  ! arbitrary large odd
    integer(int64), parameter :: mask31       = 2147483647_int64            ! 0x7FFFFFFF
    integer(int64) :: sm_state, key, w
    integer, allocatable :: state(:)
    integer :: n, i, sf

    sf = 0
    if (present(second_factor)) sf = second_factor

    if (master_seed_loaded) then
        key = ieor(master_seed_u64, int(thread_num, int64) * gamma)
        key = ieor(key,             int(sf,         int64) * second_salt)
    else
        key = int(20180815, int64) + int(thread_num, int64) + int(sf, int64) * gamma
    endif

    sm_state = key
    call splitmix64_next(sm_state, w)   ! initial diffusion of the key

    call random_seed(size=n)
    allocate(state(n))
    do i = 1, n
        call splitmix64_next(sm_state, w)
        state(i) = int(iand(w, mask31), kind(state(i)))
    enddo
    call random_seed(put=state)
    deallocate(state)

end subroutine initialise_seed_for_random_numbers_for_each_thread

end module utils_parallelised_particle_pushing_mod
