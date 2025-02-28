module utils_parallelised_particle_pushing_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine print_progress(num_particles,kpart,n)

    integer :: num_particles, kpart, n
    logical :: print_progress_for_very_particle = .false.

    if (print_progress_for_very_particle) then
        print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
    else
        if (modulo(kpart,int(num_particles/10)).eq.0) then
            print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
        endif
    endif    

end subroutine print_progress

subroutine handle_lost_particles(local_counter, boole_particle_lost)

    use boltzmann_types_mod, only: counter_t

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

subroutine initialise_loop_variables(l, n, v0, local_counter,particle_status,t,local_tetr_moments,x,vpar,vperp)

use boltzmann_types_mod, only: in, counter_t, start, time_t, particle_status_t
use constants, only: ev2erg
use tetra_physics_mod, only: particle_mass

integer, intent(in) :: l, n
type(counter_t) :: local_counter
type(particle_status_t) :: particle_status
type(time_t) :: t
logical :: boole_particle_lost
real(dp) :: t_step, t_confined, pitchpar, v, vpar, vperp, v0
real(dp), dimension(3) :: x
complex(dp), dimension(:,:) :: local_tetr_moments


call set_counter_zero(local_counter)
particle_status%lost = .false.
particle_status%initialized = .false.
particle_status%exit = .false.
t%step = in%time_step
t%confined = 0.0_dp
if (l.eq.1) local_tetr_moments = 0.0_dp
pitchpar = start%pitch(n)
x = start%x(:,n)
vpar = pitchpar * v0
vperp = sqrt(v0**2-vpar**2)
if (in%boole_boltzmann_energies) then
v = sqrt(start%energy(n)*ev2erg*2/particle_mass)
vpar = pitchpar * v
vperp = sqrt(v**2-vpar**2)
endif

end subroutine initialise_loop_variables

subroutine add_local_tetr_moments_to_output(local_tetr_moments)

    use boltzmann_types_mod, only: moment_specs, output
    use tetra_grid_mod, only: ntetr
    
    complex(dp), dimension(:,:), intent(in) :: local_tetr_moments
    integer :: k, n_prisms
    
    n_prisms = ntetr/3
    
    output%tetr_moments = output%tetr_moments + local_tetr_moments
    output%prism_moments = output%prism_moments + &
        & (local_tetr_moments(:,(/(1+3*k,k = 0,n_prisms-1)/)) + &
        &  local_tetr_moments(:,(/(2+3*k,k = 0,n_prisms-1)/)) + &
        &  local_tetr_moments(:,(/(3+3*k,k = 0,n_prisms-1)/)))
    if (moment_specs%boole_squared_moments) then
        output%prism_moments_squared = output%prism_moments_squared + &
            & (local_tetr_moments(:,(/(1+3*k,k = 0,n_prisms-1)/)) + &
            &  local_tetr_moments(:,(/(2+3*k,k = 0,n_prisms-1)/)) + &
            &  local_tetr_moments(:,(/(3+3*k,k = 0,n_prisms-1)/)))**2
    endif
    
end subroutine add_local_tetr_moments_to_output

subroutine set_counter_zero(counter)

    use boltzmann_types_mod, only: counter_t
    
    type(counter_t), intent(inout) :: counter
    
    counter%lost_particles = 0
    counter%lost_inside = 0
    counter%tetr_pushings = 0
    counter%phi_0_mappings = 0
    
end subroutine set_counter_zero

subroutine add_local_counter_to_counter(local_counter)

    use boltzmann_types_mod, only: counter_t, counter
    
    type(counter_t), intent(in) :: local_counter
    
    counter%lost_particles = counter%lost_particles + local_counter%lost_particles
    counter%lost_inside = counter%lost_inside + local_counter%lost_inside
    counter%tetr_pushings = counter%tetr_pushings + local_counter%tetr_pushings
    counter%phi_0_mappings = counter%phi_0_mappings + local_counter%phi_0_mappings
    
end subroutine add_local_counter_to_counter

subroutine carry_out_collisions(i, n, v0, t, x, vpar, vperp, ind_tetr, iface)

    use boltzmann_types_mod, only: in, c, time_t, start
    use collis_ions, only: stost
    use find_tetra_mod, only: find_tetra
    use collis_ions, only: collis_init
    use tetra_physics_mod, only: particle_mass,particle_charge
    use constants, only: echarge,amp
    
    integer, intent(in) :: i, n
    real(dp), intent(in) :: v0
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(inout) :: vpar, vperp
    type(time_t) :: t
    integer :: ind_tetr, iface
    
    real(dp), dimension(5) :: zet
    real(dp), dimension(3) :: randnum
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat,vpar_background
    real(dp) :: m0, z0, vpar_save, vperp_save, delta_epsilon, delta_vpar, vpar_mat_save, vpar_mat
    integer :: err, j
    real(dp) :: particle_to_background_coupling_strength = 1.0_dp
    
    allocate(efcolf(c%n))
    allocate(velrat(c%n))
    allocate(enrat(c%n))
    allocate(vpar_background(c%n))
    
    if (i.eq.1) call find_tetra(x,vpar,vperp,ind_tetr,iface)
    if (.not.(ind_tetr.eq.-1)) then
        if (in%boole_preserve_energy_and_momentum_during_collisions) then
            !use temp and dens to call collis_init to compute efcolf, velrat enrat and vpar_background
            do j = 1,c%n
                !> if statement because electron density will be calculated in collis init
                if (j.lt.c%n) c%dens(j) = c%dens_mat(j,ind_tetr)
                c%temp(j) = c%temp_mat(j,ind_tetr)
            enddo
            m0 = particle_mass/amp
            z0 = particle_charge/echarge
            call collis_init(m0,z0,c%mass_num,c%charge_num,c%dens,c%temp,in%energy_eV,v0,efcolf,velrat,enrat)
            vpar_save = vpar
            vperp_save = vperp
        else
            efcolf = c%efcolf_mat(:,ind_tetr)
            velrat = c%velrat_mat(:,ind_tetr)
            enrat = c%enrat_mat(:,ind_tetr)
        endif
        vpar_background = c%vpar_mat(:,ind_tetr)

        vpar = vpar - vpar_background(1)
        !since vpar_background actually has num_background_particles entries, consider giving it as an extra
        !optional input variable to stost, before randnum (maybe also check if radnum could then be set by 
        !randnum = variable even if vpar_background is not set and other variables are not set by name indexing)
        !since it came up when writing these lines: replace expressions like
        zet(1:3) = x !spatial position
        zet(4) = sqrt(vpar**2+vperp**2)/v0 !normalized velocity module 
        zet(5) = vpar/sqrt(vpar**2+vperp**2) !pitch parameter
        
        if (in%boole_precalc_collisions) then
            randnum = c%randcol(n,mod(i-1,c%randcoli)+1,:) 
            call stost(efcolf,velrat,enrat,zet,t%step,1,err,(in%time_step-t%confined)*v0,randnum)
        else
            call stost(efcolf,velrat,enrat,zet,t%step,1,err,(in%time_step-t%confined)*v0)
        endif
        
        vpar = zet(5)*zet(4)*v0+vpar_background(1)
        vperp = sqrt(1-zet(5)**2)*zet(4)*v0
        
        !optionally still change particle_mass, particle_charge and cm_over_e, e.g.:
        !particle_charge = particle_charge + echarge
        !particle_mass = particle_mass + ame - amp
        !cm_over_e = clight*particle_mass/particle_charge

        if (in%boole_preserve_energy_and_momentum_during_collisions) then
            delta_vpar = vpar - vpar_save
            delta_epsilon = particle_mass/2*(vpar**2 + vperp**2 - vpar_save**2 - vperp_save**2)

            vpar_mat_save = c%vpar_mat(1,ind_tetr)
            
            !$omp critical
            c%vpar_mat(1,ind_tetr) = vpar_mat_save - &
                                    0.01_dp*c%weight_factor*start%weight(n)*delta_vpar*particle_to_background_coupling_strength
            vpar_mat = c%vpar_mat(1,ind_tetr)
            c%temp_mat(1,ind_tetr) = c%temp_mat(1,ind_tetr) + particle_mass/3*(vpar_mat_save**2 - vpar_mat**2) - &
                                     0.01_dp*c%weight_factor*start%weight(n)*delta_epsilon*particle_to_background_coupling_strength
            !$omp end critical
        endif

    endif

end subroutine carry_out_collisions

subroutine update_exit_data(boole_particle_lost,t_confined,x,vpar,vperp,i,n,phi_0_mappings)

    use boltzmann_types_mod, only: exit_data, in

    integer, intent(in)    :: i, n
    integer, optional      :: phi_0_mappings
    real(dp), dimension(3) :: x
    real(dp)               :: t_confined, vpar, vperp
    logical                :: boole_particle_lost

    if (boole_particle_lost.eqv..true.) exit_data%lost(n) = 1
    if (boole_particle_lost.eqv..false.) exit_data%lost(n) = 0
    exit_data%t_confined(n) = t_confined
    exit_data%x(:,n) = x
    exit_data%vpar(n) = vpar
    exit_data%vperp(n) = vperp 
    exit_data%integration_step(n) = i
    if(present(phi_0_mappings)) exit_data%phi_0_mappings(n) = phi_0_mappings

end subroutine update_exit_data

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

subroutine initialise_seed_for_random_numbers_for_each_thread(thread_num)
    !This routine sets an individual seed for random number generation in each thread. It does so by adding the thread number
    !to a given array of integers and using the sum as a put argument of random_seed. Since the seed has very low entropy (every
    !value of the array is identical), the first random numbers produced are likely to be very non-random. Thus, n random
    !numbers are generated to get rid of these potentially corrupted numnbers (compare with 
    !https://stackoverflow.com/questions/51893720/correctly-setting-random-seeds-for-repeatability, also check
    !https://stats.stackexchange.com/questions/354373/what-exactly-is-a-seed-in-a-random-number-generator)

    integer, intent(in) :: thread_num
    real(dp) :: randnum
    integer :: i,n, state(33)

    n = 1000
    state = 20180815

    call random_seed(put=state+thread_num)

    do i = 1,n
        call random_number(randnum)
    enddo

end subroutine initialise_seed_for_random_numbers_for_each_thread

end module utils_parallelised_particle_pushing_mod