module orbit_timestep_gorilla_supporting_functions_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine categorize_lost_particles(ind_tetr,x,local_counter,t_remain,t_remain_out)

    use tetra_physics_mod, only: tetra_physics
    use boltzmann_types_mod, only: counter_t, pflux

    integer, intent(in) :: ind_tetr
    real(dp), dimension(3), intent(in) :: x
    real(dp), intent(in) :: t_remain
    type(counter_t), intent(inout) :: local_counter
    real(dp), intent(out), optional              :: t_remain_out
    real(dp) :: aphi

    !print *, 'WARNING: Particle lost.'
    aphi = tetra_physics(ind_tetr)%Aphi1+sum(tetra_physics(ind_tetr)%gAphi*x)
    if (abs(aphi-pflux%max).lt.abs(aphi-pflux%min)) then
        local_counter%lost_outside = local_counter%lost_outside + 1
    endif
    if (abs(aphi-pflux%max).gt.abs(aphi-pflux%min)) then
        local_counter%lost_inside = local_counter%lost_inside + 1
    endif
    if( present(t_remain_out)) then
        t_remain_out = t_remain
    endif

end subroutine categorize_lost_particles

subroutine update_local_tetr_moments(local_tetr_moments,ind_tetr,n,optional_quantities)

    use boltzmann_types_mod, only: moment_specs, start
    use gorilla_settings_mod, only: optional_quantities_type

    type(optional_quantities_type), intent(in)   :: optional_quantities
    integer, intent(in)                          :: ind_tetr, n
    complex(dp), dimension(:,:), intent (inout)  :: local_tetr_moments
    integer                                      :: m

    do m = 1,moment_specs%n_moments
        select case(moment_specs%moments_selector(m))
            case(1)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                                & start%weight(n)*optional_quantities%t_hamiltonian
            case(2)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                                & start%weight(n)*optional_quantities%gyrophase*start%jperp(n)
            case(3)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                                & start%weight(n)*optional_quantities%vpar_int
            case(4)
                local_tetr_moments(m,ind_tetr) = local_tetr_moments(m,ind_tetr) + &
                                                                & start%weight(n)*optional_quantities%vpar2_int
        end select
    enddo

end subroutine update_local_tetr_moments

subroutine initialize_constants_of_motion(vperp,z_save,ind_tetr,perpinv)

    use pusher_tetra_rk_mod, only: initialize_const_motion_rk
    use pusher_tetra_poly_mod, only: initialize_const_motion_poly
    use gorilla_settings_mod, only: ipusher
    use supporting_functions_mod, only: bmod_func

    real(dp), intent(in)               :: vperp
    real(dp), dimension(3), intent(in) :: z_save
    integer, intent(in)                :: ind_tetr
    real(dp), intent(out)              :: perpinv

    perpinv=-0.5d0*vperp**2/bmod_func(z_save,ind_tetr)
             
    !Initialize constants of motion in particle-private module
    select case(ipusher)
        case(1)
            call initialize_const_motion_rk(perpinv,perpinv**2)
        case(2)
            call initialize_const_motion_poly(perpinv,perpinv**2)
    end select

end subroutine initialize_constants_of_motion

subroutine calc_particle_weights_and_jperp(n,z_save,vpar,vperp,ind_tetr)

    use boltzmann_types_mod, only: u, pflux, start
    use tetra_physics_mod, only: tetra_physics,particle_mass,particle_charge,cm_over_e
    use constants, only: ev2erg
    use volume_integrals_and_sqrt_g_mod, only: sqrt_g
    use supporting_functions_mod, only: bmod_func

    real(dp), intent(in) :: vpar, vperp
    real(dp), dimension(3), intent(in) :: z_save
    integer, intent(in) :: n,ind_tetr
    real(dp) :: local_poloidal_flux, phi_elec_func, temperature
    real(dp) :: r, phi, z

    r = z_save(1)
    phi = z_save(2)
    z = z_save(3)

    if (.not.u%boole_refined_sqrt_g) start%weight = start%weight*r
    if (u%boole_refined_sqrt_g) then
        start%weight(n) = start%weight(n)* (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
                                        &  (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))               
    endif

    if (u%boole_linear_density_simulation.or.u%boole_linear_temperature_simulation) then
        local_poloidal_flux = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*z_save) + &
                                        & cm_over_e*vpar*&
                                        & (tetra_physics(ind_tetr)%h2_1+sum(tetra_physics(ind_tetr)%gh2*z_save))
    endif
    if (u%boole_linear_density_simulation) then
        start%weight(n) = start%weight(n)*(pflux%max*1.1-local_poloidal_flux)/(pflux%max*1.1)
    endif

    if (u%boole_boltzmann_energies) then
        !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts have been added before)
        phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z_save)
        if (.not. u%boole_linear_temperature_simulation) then
            start%weight(n) =start%weight(n)*sqrt(start%energy(n)*ev2erg)/(u%energy_eV*ev2erg)**1.5* &
                        & exp(-(start%energy(n)*ev2erg+particle_charge*phi_elec_func)/(u%energy_eV*ev2erg))
        else
            temperature = u%energy_eV*ev2erg*(pflux%max*1.1-local_poloidal_flux)/(pflux%max*1.1)
            start%weight(n) = start%weight(n)*sqrt(start%energy(n)*ev2erg)/temperature**1.5* &
            & exp(-(start%energy(n)*ev2erg+particle_charge*phi_elec_func)/temperature)
        endif
    endif

    start%jperp(n) = particle_mass*vperp**2*cm_over_e/(2*bmod_func(z_save,ind_tetr))*(-1) !-1 because of negative gyrophase

end subroutine calc_particle_weights_and_jperp

end module orbit_timestep_gorilla_supporting_functions_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module main_routine_supporting_functions_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine print_progress(num_particles,kpart,n)

    integer :: num_particles, kpart, n

    ! if (num_particles.gt.10) then
    !     if (modulo(kpart,int(num_particles/10)).eq.0) then
    !         print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
    !     endif
    ! else 
    !     print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
    ! endif

    print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()

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

subroutine initialise_loop_variables(l, n, v0, local_counter,boole,t,local_tetr_moments,x,vpar,vperp)

use boltzmann_types_mod, only: u, counter_t, start, time_t, boole_t
use constants, only: ev2erg
use tetra_physics_mod, only: particle_mass

integer, intent(in) :: l, n
type(counter_t) :: local_counter
type(boole_t) :: boole
type(time_t) :: t
logical :: boole_particle_lost
real(dp) :: t_step, t_confined, pitchpar, v, vpar, vperp, v0
real(dp), dimension(3) :: x
complex(dp), dimension(:,:) :: local_tetr_moments


call set_local_counter_zero(local_counter)
boole%lost = .false.
boole%initialized = .false.
boole%exit = .false.
t%step = u%time_step
t%confined = 0
if (l.eq.1) local_tetr_moments = 0
pitchpar = start%pitch(n)
x = start%x(:,n)
vpar = pitchpar * v0
vperp = sqrt(v0**2-vpar**2)
if (u%boole_boltzmann_energies) then
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

subroutine normalise_prism_moments_and_prism_moments_squared

    use boltzmann_types_mod, only: moment_specs, output, u
    
    integer :: n
    
    
    do n = 1,moment_specs%n_moments
    output%prism_moments(n,:) = output%prism_moments(n,:)/(output%prism_volumes*u%time_step*u%n_particles)
    if (moment_specs%boole_squared_moments) then
    output%prism_moments_squared(n,:) = output%prism_moments_squared(n,:)/ &
                   (output%prism_volumes**2*u%time_step**2*u%n_particles)
    endif
    if (u%boole_refined_sqrt_g) then
    output%prism_moments(n,:) = output%prism_moments(n,:)*output%prism_volumes/output%refined_prism_volumes
    if (moment_specs%boole_squared_moments) then
    output%prism_moments_squared(n,:) = output%prism_moments_squared(n,:)* &
           output%prism_volumes**2/output%refined_prism_volumes**2
    endif
    endif
    enddo
    
end subroutine normalise_prism_moments_and_prism_moments_squared

subroutine set_moment_specifications

    use gorilla_settings_mod, only: boole_array_optional_quantities
    use tetra_grid_settings_mod, only: grid_size
    use tetra_grid_mod, only: ntetr
    use boltzmann_types_mod, only: moment_specs, u
    
    integer :: i, n_prisms
    
    n_prisms = ntetr/3
    
    moment_specs%boole_squared_moments = u%boole_squared_moments
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

subroutine set_local_counter_zero(counter)

    use boltzmann_types_mod, only: counter_t
    
    type(counter_t), intent(inout) :: counter
    
    counter%lost_particles = 0
    counter%lost_inside = 0
    counter%lost_outside = 0
    counter%tetr_pushings = 0
    counter%phi_0_mappings = 0
    
end subroutine set_local_counter_zero

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

subroutine add_local_counter_to_counter(local_counter)

    use boltzmann_types_mod, only: counter_t, counter
    
    type(counter_t), intent(in) :: local_counter
    
    counter%lost_particles = counter%lost_particles + local_counter%lost_particles
    counter%lost_inside = counter%lost_inside + local_counter%lost_inside
    counter%lost_outside = counter%lost_outside + local_counter%lost_outside
    counter%tetr_pushings = counter%tetr_pushings + local_counter%tetr_pushings
    counter%phi_0_mappings = counter%phi_0_mappings + local_counter%phi_0_mappings
    
end subroutine add_local_counter_to_counter

subroutine get_ipert()

    use field_mod, only: ipert
    
    integer :: ipert_unit
    
    open(newunit = ipert_unit, file='field_divB0.inp')
    read(ipert_unit,*) ipert        ! 0=eq only, 1=vac, 2=vac+plas no derivatives,
    close(ipert_unit)
    
end subroutine get_ipert

subroutine calc_poloidal_flux(verts)

    use boltzmann_types_mod, only: pflux
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: ntetr, tetra_grid
    
    real(dp), dimension(:,:), intent(in) :: verts
    integer :: i
    
    pflux%max = 0
    pflux%min = tetra_physics(1)%Aphi1
    do i = 1, ntetr
    pflux%max = max(pflux%max,tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
    & (verts([1,2,3],tetra_grid(i)%ind_knot(4))-verts([1,2,3],tetra_grid(i)%ind_knot(1)))))
    pflux%min = min(pflux%min,tetra_physics(i)%Aphi1)
    enddo
    
end subroutine calc_poloidal_flux

subroutine calc_starting_conditions(verts)

    use boltzmann_types_mod, only: u
    
    real(dp), dimension(:,:), allocatable, intent(out)     :: verts
    real(dp), dimension(:,:), allocatable                  :: rand_matrix

    call set_verts_and_coordinate_limits(verts)
    allocate(rand_matrix(5,u%num_particles))
    call RANDOM_NUMBER(rand_matrix)
    call allocate_start_type
    call set_starting_positions(rand_matrix)
    call set_rest_of_start_type(rand_matrix)

end subroutine calc_starting_conditions

subroutine set_starting_positions(rand_matrix)

    use boltzmann_types_mod, only: u, start, g
    use tetra_physics_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind
    use constants, only: pi

    real(dp), dimension(:,:), intent(in) :: rand_matrix

    !compute starting conditions
    if (u%boole_point_source) then
        if (grid_kind.eq.2) then
            start%x(1,:) = 160
            start%x(2,:) = 0
            start%x(3,:) = 70
        elseif (grid_kind.eq.4) then
            start%x(1,:) = 205
            start%x(2,:) = 0
            start%x(3,:) = 0
        endif
        if (coord_system.eq.2) print*, 'error: point source is only implemented for cylindrical coordinate system'
    else
        start%x(g%ind_a,:) = g%amin + (g%amax - g%amin)*rand_matrix(g%ind_a,:) !r in cylindrical, s in flux coordinates
        start%x(g%ind_b,:) = 2*pi*rand_matrix(g%ind_b,:) !phi in cylindrical and flux coordinates
        start%x(g%ind_c,:) = g%cmin + (g%cmax - g%cmin)*rand_matrix(g%ind_c,:) !z in cylindrical, theta in flux coordinates
    endif

end subroutine set_starting_positions

subroutine set_rest_of_start_type(rand_matrix)

    use boltzmann_types_mod, only: u, start, g
    use constants, only: pi, ev2erg

    real(dp), dimension(:,:), intent(in) :: rand_matrix

    start%pitch(:) = 2*rand_matrix(4,:)-1 !pitch parameter
    start%energy = u%energy_eV
    start%weight = u%density*(g%amax-g%amin)*(g%cmax-g%cmin)*2*pi
    if (u%boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start%energy = 5*u%energy_eV*rand_matrix(5,:) !boltzmann energy distribution
        start%weight =  start%weight*10/sqrt(pi)*u%energy_eV*ev2erg
    endif
    
    if (u%boole_antithetic_variate) then
        start%x(:,1:u%num_particles:2) = start%x(:,2:u%num_particles:2)
        start%pitch(1:u%num_particles:2) = -start%pitch(2:u%num_particles:2)
        start%energy(1:u%num_particles:2) = start%energy(2:u%num_particles:2)
    endif

end subroutine set_rest_of_start_type

subroutine set_verts_and_coordinate_limits(verts)

    use tetra_physics_mod, only: coord_system
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, nvert
    use boltzmann_types_mod, only: g

    real(dp), dimension(:,:), allocatable, intent(out)     :: verts

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

end subroutine set_verts_and_coordinate_limits

subroutine allocate_start_type

    use boltzmann_types_mod, only: start, u

    allocate(start%x(3,u%num_particles))
    allocate(start%pitch(u%num_particles))
    allocate(start%energy(u%num_particles))
    allocate(start%weight(u%num_particles))
    allocate(start%jperp(u%num_particles))

end subroutine allocate_start_type

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

subroutine calc_collision_coefficients_for_all_tetrahedra(v0)

    use boltzmann_types_mod, only: u, c
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
    
    do i = 1, ntetr
    do j = 1,c%n
    !> if statement because electron density will be calculated in collis init
    if (j.lt.c%n) c%dens(j) = c%dens_mat(j,i)
    c%temp(j) = c%temp_mat(j,i)
    enddo
    call collis_init(m0,z0,c%mass_num,c%charge_num,c%dens,c%temp,u%energy_eV,v0,efcolf,velrat,enrat)
    c%efcolf_mat(:,i) = efcolf
    c%velrat_mat(:,i) = velrat
    c%enrat_mat(:,i) = enrat
    enddo
    
    if (u%boole_precalc_collisions) then
    allocate(c%randcol(u%num_particles,c%randcoli,3))
    call RANDOM_NUMBER(c%randcol)
    !3.464102 = sqrt(12), this creates a random number with zero average and unit variance
    c%randcol(:,:,1:2:3) =  3.464102*(c%randcol(:,:,1:2:3)-.5)
    endif
end subroutine calc_collision_coefficients_for_all_tetrahedra

subroutine carry_out_collisions(i, n, v0, t, x, vpar, vperp, ind_tetr, iface)

    use boltzmann_types_mod, only: u, c, time_t
    use collis_ions, only: stost
    use find_tetra_mod, only: find_tetra
    
    integer, intent(in) :: i, n
    real(dp), intent(in) :: v0
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(inout) :: vpar, vperp
    type(time_t) :: t
    integer :: ind_tetr, iface
    
    real(dp), dimension(5) :: zet
    real(dp), dimension(3) :: randnum
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat,vpar_background
    integer :: err
    
    allocate(efcolf(c%n))
    allocate(velrat(c%n))
    allocate(enrat(c%n))
    allocate(vpar_background(c%n))
    
    if (i.eq.1) call find_tetra(x,vpar,vperp,ind_tetr,iface)
    if (.not.(ind_tetr.eq.-1)) then
    efcolf = c%efcolf_mat(:,ind_tetr)
    velrat = c%velrat_mat(:,ind_tetr)
    enrat = c%enrat_mat(:,ind_tetr)
    vpar_background = c%vpar_mat(:,ind_tetr)
    
    vpar = vpar - vpar_background(1)
    !since vpar_background actually has num_background_particles entries, consider giving it as an extra
    !optional input variable to stost, before randnum (maybe also check if radnum could then be set by 
    !randnum = variable eve if vpar_background is not set and other variables are not set by name indexing)
    !since it came up when writing these lines: replace expressions like
    zet(1:3) = x !spatial position
    zet(4) = sqrt(vpar**2+vperp**2)/v0 !normalized velocity module 
    zet(5) = vpar/sqrt(vpar**2+vperp**2) !pitch parameter
    
    if (u%boole_precalc_collisions) then
    randnum = c%randcol(n,mod(i-1,c%randcoli)+1,:) 
    call stost(efcolf,velrat,enrat,zet,t%step,1,err,(u%time_step-t%confined)*v0,randnum)
    else
    call stost(efcolf,velrat,enrat,zet,t%step,1,err,(u%time_step-t%confined)*v0)
    endif
    
    vpar = zet(5)*zet(4)*v0+vpar_background(1)
    vperp = sqrt(1-zet(5)**2)*zet(4)*v0
    
    !optionally still change particle_mass, particle_charge and cm_over_e, e.g.:
    !particle_charge = particle_charge + echarge
    !particle_mass = particle_mass + ame - amp 
    !cm_over_e = clight*particle_mass/particle_charge
    endif
    
end subroutine carry_out_collisions

subroutine initialize_exit_data

    use boltzmann_types_mod, only: u, exit_data

    allocate(exit_data%lost(u%num_particles))
    allocate(exit_data%t_confined(u%num_particles))
    allocate(exit_data%x(3,u%num_particles))
    allocate(exit_data%vpar(u%num_particles))
    allocate(exit_data%vperp(u%num_particles))
    allocate(exit_data%phi_0_mappings(u%num_particles))
    allocate(exit_data%integration_step(u%num_particles))

    exit_data%lost = 0
    exit_data%t_confined = 0.0d0
    exit_data%x = 0.0d0
    exit_data%vpar = 0.0d0 
    exit_data%vperp = 0.0d0 
    exit_data%integration_step = 0
    exit_data%phi_0_mappings = 0


end subroutine initialize_exit_data

subroutine update_exit_data(boole_particle_lost,t_confined,x,vpar,vperp,i,n,phi_0_mappings)

    use boltzmann_types_mod, only: exit_data, u

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

end module main_routine_supporting_functions_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module boltzmann_writing_data_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    public :: write_data_to_files, name_files, unlink_files

contains

subroutine write_data_to_files

    use boltzmann_types_mod, only: filenames, u, output, moment_specs

    if (u%boole_write_vertex_indices) call write_vertex_indices

    if (u%boole_write_vertex_coordinates) call write_vertex_coordinates

    if (u%boole_write_prism_volumes) call write_prism_volumes

    if (u%boole_write_refined_prism_volumes) call write_refined_prism_volumes

    if (u%boole_write_boltzmann_density) call write_boltzmann_densities

    if (u%boole_write_electric_potential) call write_electric_potential

    if (u%boole_write_moments) then
        if (moment_specs%n_moments.gt.0) then
            call write_moments
        else
            print*, "Error: moments are not written to file because no moment was computed. Turn computation of moments on in &
                     gorilla.inp."
        endif
    endif

    if (u%boole_write_fourier_moments) then
        if (moment_specs%n_moments.gt.0) then
            call write_fourier_moments
        else
            print*, "Error: Fourier moments are not written to file because no moment was computed (and thus also no fourier &
                     moment). Turn computation of moments on in gorilla.inp."
        endif

    endif

    if (u%boole_write_exit_data) call write_exit_data

end subroutine write_data_to_files

subroutine write_vertex_indices

    use tetra_grid_mod, only: ntetr, tetra_grid
    use boltzmann_types_mod, only: filenames

    integer :: vi_unit
    integer :: i

    open(newunit = vi_unit, file = filenames%vertex_indices)
    do i=1, ntetr
        write(vi_unit, *) tetra_grid(i)%ind_knot([1, 2, 3, 4])
    end do
    close(vi_unit)

end subroutine write_vertex_indices

subroutine write_vertex_coordinates

    use tetra_physics_mod, only: coord_system
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, nvert
    use boltzmann_types_mod, only: filenames

    integer :: vc_unit
    integer :: i

    open(newunit = vc_unit, file = filenames%vertex_coordinates)
    101 format(1000(e21.14,x))
    if (coord_system.eq.1) then
        ![R,phi,Z]: Write vertex coordinates to file
        do i=1, nvert
            write(vc_unit,101) verts_rphiz(1, i), verts_rphiz(2, i), verts_rphiz(3, i)
        end do
    elseif (coord_system.eq.2) then
        ![s,theta,phi]: Write vertex coordinates to file
        do i=1, nvert
            write(vc_unit,101) verts_sthetaphi(1, i), verts_sthetaphi(2, i), verts_sthetaphi(3, i)
        end do
    endif
    close(vc_unit)

end subroutine write_vertex_coordinates

subroutine write_prism_volumes

    use boltzmann_types_mod, only: filenames, output

    integer :: pv_unit

    open(newunit = pv_unit, file = filenames%prism_volumes)
    write(pv_unit,'(ES20.10E4)') output%prism_volumes
    close(pv_unit)

end subroutine write_prism_volumes

subroutine write_refined_prism_volumes

    use boltzmann_types_mod, only: filenames, output

    integer :: rpv_unit

    open(newunit = rpv_unit, file = filenames%refined_prism_volumes)
    write(rpv_unit,'(ES20.10E4)') output%refined_prism_volumes
    close(rpv_unit)

end subroutine write_refined_prism_volumes

subroutine write_boltzmann_densities

    use boltzmann_types_mod, only: filenames, output

    integer :: bd_unit

    open(newunit = bd_unit, file = filenames%boltzmann_density)
    write(bd_unit,'(ES20.10E4)') output%boltzmann_density
    close(bd_unit)

end subroutine write_boltzmann_densities

subroutine write_electric_potential

    use boltzmann_types_mod, only: filenames, output

    integer :: epv_unit

    open(newunit = epv_unit, file = filenames%electric_potential)
    write(epv_unit,'(ES20.10E4)') output%electric_potential
    close(epv_unit)

end subroutine write_electric_potential

subroutine write_moments

    use tetra_grid_mod, only: ntetr
    use boltzmann_types_mod, only: moment_specs, filenames, output

    integer :: p_moments_unit, pmss_unit, t_moments_unit
    integer :: l, i
    integer :: n_prisms

    open(newunit = p_moments_unit, file = filenames%prism_moments)
    open(newunit = pmss_unit, file = filenames%prism_moments_summed_squares)
    open(newunit = t_moments_unit, file = filenames%tetr_moments)

    n_prisms = ntetr/3

    if (moment_specs%n_moments.gt.0) then
        do l = 1,n_prisms
            do i = 1,moment_specs%n_moments - 1
                write(p_moments_unit,'(2ES20.10E4)',advance="no") real(output%prism_moments(i,l)), aimag(output%prism_moments(i,l))
            enddo
                write(p_moments_unit,'(2ES20.10E4)') real(output%prism_moments(moment_specs%n_moments,l)), &
                                                     aimag(output%prism_moments(moment_specs%n_moments,l))
        enddo
        if (moment_specs%boole_squared_moments) then
            do l = 1,n_prisms
                do i = 1,moment_specs%n_moments - 1
                    write(pmss_unit,'(2ES20.10E4)',advance="no") real(output%prism_moments_squared(i,l)), &
                                                                    & aimag(output%prism_moments_squared(i,l))
                enddo
                    write(pmss_unit,'(2ES20.10E4)') real(output%prism_moments_squared(moment_specs%n_moments,l)), &
                                                    aimag(output%prism_moments_squared(moment_specs%n_moments,l))
            enddo
        endif
        do l = 1,ntetr
            do i = 1,moment_specs%n_moments - 1
                write(t_moments_unit,'(2ES20.10E4)',advance="no") real(output%tetr_moments(i,l)), aimag(output%tetr_moments(i,l))
            enddo
                write(t_moments_unit,'(2ES20.10E4)') real(output%tetr_moments(moment_specs%n_moments,l)), &
                                                     aimag(output%tetr_moments(moment_specs%n_moments,l))
        enddo
    endif

    close(p_moments_unit)
    close(pmss_unit)
    close(t_moments_unit)

end subroutine write_moments

subroutine write_fourier_moments

    use tetra_grid_settings_mod, only: grid_size
    use boltzmann_types_mod, only: moment_specs, output, filenames

    integer :: fm_unit, l, i

    open(newunit = fm_unit, file = filenames%fourier_moments)
    do l = 1,moment_specs%n_triangles
        do i = 1,moment_specs%n_fourier_modes-1
            write(fm_unit,'(2ES20.10E4)',advance="no") real(output%moments_in_frequency_space(1,l,i)), &
                                                       aimag(output%moments_in_frequency_space(1,l,i))
        enddo
            write(fm_unit,'(2ES20.10E4)') real(output%moments_in_frequency_space(1,l,moment_specs%n_fourier_modes)), &
                                          aimag(output%moments_in_frequency_space(1,l,moment_specs%n_fourier_modes))
    enddo
    close(fm_unit)

end subroutine write_fourier_moments

subroutine write_exit_data

    use boltzmann_types_mod, only: exit_data, filenames, u

    integer :: ed_unit, i

    open(newunit = ed_unit, file = filenames%exit_data)
    do i = 1,u%num_particles
        write(ed_unit,*) i, exit_data%lost(i), exit_data%t_confined(i), exit_data%x(:,i), exit_data%vpar(i), &
                        exit_data%vperp(i), exit_data%integration_step(i), exit_data%phi_0_mappings(i)
    enddo
    close (ed_unit)

end subroutine write_exit_data

subroutine name_files

    use boltzmann_types_mod, only: filenames

    filenames%poincare_maps = 'poincare_maps.dat'
    filenames%prism_moments = 'prism_moments.dat'
    filenames%prism_moments_summed_squares = 'prism_moments_summed_squares.dat'
    filenames%vertex_coordinates = 'vertex_coordinates.dat'
    filenames%vertex_indices = 'vertex_indices.dat'
    filenames%prism_volumes = 'prism_volumes.dat'
    filenames%fourier_moments = 'fourier_moments.dat'
    filenames%refined_prism_volumes = 'refined_prism_volumes.dat'
    filenames%electric_potential = 'electric_potential.dat'
    filenames%boltzmann_density = 'boltzmann_density.dat'
    filenames%divertor_intersections = 'divertor_intersections.dat'
    filenames%tetr_moments = 'tetr_moments.dat'
    filenames%exit_data = 'exit_data.dat'

end subroutine name_files

subroutine unlink_files

    use boltzmann_types_mod, only: filenames

    call unlink(filenames%poincare_maps)
    call unlink(filenames%prism_moments)
    call unlink(filenames%prism_moments_summed_squares)
    call unlink(filenames%vertex_coordinates)
    call unlink(filenames%vertex_indices)
    call unlink(filenames%prism_volumes)
    call unlink(filenames%fourier_moments)
    call unlink(filenames%refined_prism_volumes)
    call unlink(filenames%electric_potential)
    call unlink(filenames%boltzmann_density)
    call unlink(filenames%divertor_intersections)
    call unlink(filenames%tetr_moments)
    call unlink(filenames%exit_data)

end subroutine unlink_files

end module boltzmann_writing_data_mod