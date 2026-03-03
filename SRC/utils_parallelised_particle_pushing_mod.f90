module utils_parallelised_particle_pushing_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine print_progress(num_particles,kpart,n)

    integer :: num_particles, kpart, n
    logical :: print_progress_for_every_particle = .false.

    if ((.not.print_progress_for_every_particle).and.(num_particles.gt.10)) then
        if (modulo(kpart,int(num_particles/10)).eq.0) then
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
use gorilla_applets_settings_mod, only: i_option

integer, intent(in) :: l, n
integer, intent(in), optional :: species_in
type(counter_t) :: local_counter
type(particle_status_t) :: particle_status
type(time_t) :: t
real(dp) :: t_step, t_confined, pitchpar, v, vpar, vperp, epsilon
real(dp), dimension(3) :: x
complex(dp), dimension(:,:) :: local_tetr_moments
integer :: species = 1

if (present(species_in)) species = species_in

call set_counter_zero(local_counter)
particle_status%lost = .false.
particle_status%initialized = .false.
particle_status%exit = .false.
t%step = in%time_step
if (present(species_in)) t%step = start%t(species)
t%confined = 0.0_dp
if (l.eq.1) local_tetr_moments = 0.0_dp
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

subroutine add_local_tetr_moments_to_output(local_tetr_moments,species_in)

    use gorilla_applets_types_mod, only: moment_specs, output
    use tetra_grid_mod, only: ntetr
    
    complex(dp), dimension(:,:), intent(in) :: local_tetr_moments
    integer, intent(in), optional :: species_in
    integer :: species = 1
    integer :: k, n_prisms

    if (present(species_in)) species = species_in
    
    n_prisms = ntetr/3
    
    output%tetr_moments(:,:,species) = output%tetr_moments(:,:,species) + local_tetr_moments
    output%prism_moments(:,:,species) = output%prism_moments(:,:,species) + &
        & (local_tetr_moments(:,(/(1+3*k,k = 0,n_prisms-1)/)) + &
        &  local_tetr_moments(:,(/(2+3*k,k = 0,n_prisms-1)/)) + &
        &  local_tetr_moments(:,(/(3+3*k,k = 0,n_prisms-1)/)))
    if (moment_specs%boole_squared_moments) then
        output%prism_moments_squared(:,:,species) = output%prism_moments_squared(:,:,species) + &
            & (local_tetr_moments(:,(/(1+3*k,k = 0,n_prisms-1)/)) + &
            &  local_tetr_moments(:,(/(2+3*k,k = 0,n_prisms-1)/)) + &
            &  local_tetr_moments(:,(/(3+3*k,k = 0,n_prisms-1)/)))**2
    endif
    
end subroutine add_local_tetr_moments_to_output

subroutine set_counter_zero(counter)

    use gorilla_applets_types_mod, only: counter_t
    
    type(counter_t), intent(inout) :: counter
    
    counter%lost_particles = 0
    counter%lost_inside = 0
    counter%tetr_pushings = 0
    counter%phi_0_mappings = 0
    counter%integration_steps = 0
    
end subroutine set_counter_zero

subroutine add_local_counter_to_counter(local_counter)

    use gorilla_applets_types_mod, only: counter_t, counter
    
    type(counter_t), intent(in) :: local_counter
    
    counter%lost_particles = counter%lost_particles + local_counter%lost_particles
    counter%lost_inside = counter%lost_inside + local_counter%lost_inside
    counter%tetr_pushings = counter%tetr_pushings + local_counter%tetr_pushings
    counter%phi_0_mappings = counter%phi_0_mappings + local_counter%phi_0_mappings
    
end subroutine add_local_counter_to_counter

subroutine carry_out_collisions(i, n, t, x, vpar, vperp, ind_tetr, iface, species_in)

    use gorilla_applets_types_mod, only: in, time_t
    use find_tetra_mod, only: find_tetra
    
    integer, intent(in) :: i, n
    integer, intent(in), optional :: species_in
    integer :: species = 1
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(inout) :: vpar, vperp
    type(time_t) :: t
    integer :: ind_tetr, iface

    if (present(species_in)) species = species_in

    if (i.eq.1) call find_tetra(x,vpar,vperp,ind_tetr,iface)
    if (.not.(ind_tetr.eq.-1)) then
        if (in%boole_preserve_energy_and_momentum_during_collisions) then
            call collisions_with_background_updates(i, n, t, x, vpar, vperp, ind_tetr, species)
        else
            call collisions_without_background_updates(i, n, t, x, vpar, vperp, ind_tetr, species)
        endif
    endif

end subroutine carry_out_collisions

subroutine collisions_with_background_updates(i, n, t, x, vpar, vperp, ind_tetr, species)

    use gorilla_applets_types_mod, only: in, c, time_t, start
    use collis_ions, only: stost
    use collis_ions, only: collis_init
    use tetra_physics_mod, only: particle_mass,particle_charge
    use constants, only: echarge,amp
    use gorilla_applets_settings_mod, only: i_option
    
    integer, intent(in) :: i, n, species
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(inout) :: vpar, vperp
    type(time_t) :: t
    integer :: ind_tetr
    
    real(dp), dimension(5) :: zet
    real(dp), dimension(3) :: randnum
    real(dp), dimension(1) :: m, z, dens, temp, efcolf,velrat,enrat
    real(dp) :: vpar_background
    real(dp) :: m0, z0, vpar_save, vperp_save, delta_epsilon, delta_vpar, vpar_mat_save, vpar_mat
    integer :: err, j, p, iswmode
    real(dp) ::  w_v, w_t, particle_to_background_coupling_strength, t_max

    iswmode = 1
    !1 - full operator (pitch-angle and energy scattering and drag)
    !2 - energy scattering and drag only
    !3 - drag only
    !4 - pitch-angle scattering only

    w_v = 1.0_dp
    w_t = 1.0_dp
    particle_to_background_coupling_strength = 0.0001_dp

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

        vpar_mat_save = c%vpar_mat(j,ind_tetr)
        
        !$omp critical
        c%vpar_mat(j,ind_tetr) = vpar_mat_save - &
                            c%weight_factor*start%weight(n,species)*delta_vpar/w_v*particle_to_background_coupling_strength
        vpar_mat = c%vpar_mat(j,ind_tetr)
        c%temp_mat(j,ind_tetr) = c%temp_mat(j,ind_tetr) + particle_mass/3*(vpar_mat_save**2 - vpar_mat**2) - &
                            2.0_dp/3.0_dp*c%weight_factor*start%weight(n,species)*delta_epsilon/w_t &
                            *particle_to_background_coupling_strength
        !$omp end critical
    enddo

end subroutine collisions_with_background_updates

subroutine collisions_without_background_updates(i, n, t, x, vpar, vperp, ind_tetr, species)

    use gorilla_applets_types_mod, only: in, c, time_t, start
    use collis_ions, only: stost
    use gorilla_applets_settings_mod, only: i_option
    
    integer, intent(in) :: i, n, species
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(inout) :: vpar, vperp
    type(time_t) :: t
    integer :: ind_tetr
    
    real(dp), dimension(5) :: zet
    real(dp), dimension(3) :: randnum
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat,vpar_background
    integer :: err, iswmode
    real(dp) :: t_max

    iswmode = 4
    !1 - full operator (pitch-angle and energy scattering and drag)
    !2 - energy scattering and drag only
    !3 - drag only
    !4 - pitch-angle scattering only

    allocate(efcolf(c%n))
    allocate(velrat(c%n))
    allocate(enrat(c%n))
    allocate(vpar_background(c%n))

    efcolf = c%efcolf_mat(:,ind_tetr)
    velrat = c%velrat_mat(:,ind_tetr)
    enrat = c%enrat_mat(:,ind_tetr)

    vpar_background = c%vpar_mat(:,ind_tetr)

    vpar = vpar - vpar_background(1)
    !since vpar_background actually has num_background_particles entries, consider giving it as an extra
    !optional input variable to stost, before randnum
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

    vpar = zet(5)*zet(4)*start%v0(species)+vpar_background(1)
    vperp = sqrt(1-zet(5)**2)*zet(4)*start%v0(species)

    !optionally still change particle_mass, particle_charge and cm_over_e, e.g.:
    !particle_charge = particle_charge + echarge
    !particle_mass = particle_mass + ame - amp
    !cm_over_e = clight*particle_mass/particle_charge

end subroutine collisions_without_background_updates

subroutine update_exit_data(boole_particle_lost,t_confined,x,vpar,vperp,i,n,phi_0_mappings,species_in,ind_tetr)

    use gorilla_applets_types_mod, only: exit_data, in, flux
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind

    integer, intent(in)                 :: i, n
    integer, intent(in), optional       :: phi_0_mappings, species_in, ind_tetr
    real(dp), dimension(3), intent(in)  :: x
    real(dp), intent(in)                :: t_confined, vpar, vperp
    logical, intent(in)                 :: boole_particle_lost
    integer                             :: species = 1
    real(dp), dimension(3)              :: z_local
    real(dp)                            :: local_poloidal_flux

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
            exit_data%flux_surface(n,species) = (local_poloidal_flux - flux%poloidal_min) / &
                                                 (flux%poloidal_max - flux%poloidal_min)
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
    !This routine sets an individual seed for random number generation in each thread. It does so by adding the thread number
    !to a given array of integers and using the sum as a put argument of random_seed. Since the seed has rather low entropy (values
    !of the array are multiples of each other), the first random numbers produced are likely to be non-random. Thus, n random
    !numbers are generated to get rid of these potentially corrupted numnbers (compare with 
    !https://stackoverflow.com/questions/51893720/correctly-setting-random-seeds-for-repeatability, also check
    !https://stats.stackexchange.com/questions/354373/what-exactly-is-a-seed-in-a-random-number-generator)

    integer, intent(in) :: thread_num
    integer, intent(in), optional :: second_factor
    real(dp) :: randnum
    integer :: i,n, state(33)

    n = 1000
    state = 20180815

    do i = 1, size(state)
        if (present(second_factor)) then 
            state(i) = (state(i)+thread_num+second_factor)*i
        else
            state(i) = (state(i)+thread_num)*i
        endif
    enddo

    call random_seed(put=state+thread_num)

    do i = 1,n
        call random_number(randnum)
    enddo

end subroutine initialise_seed_for_random_numbers_for_each_thread

end module utils_parallelised_particle_pushing_mod
