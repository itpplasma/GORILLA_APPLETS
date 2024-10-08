module orbit_timestep_gorilla_supporting_functions_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
   
contains

subroutine calc_and_write_poincare_mappings_and_divertor_intersections(x_save,x,n,iper_phi,local_counter,iunits,boole_t_finished)

    use boltzmann_types_mod, only: counter_t, iunits_t

    real(dp), dimension(3), intent(in)  :: x, x_save
    integer, intent(in)                 :: n, iper_phi
    type(iunits_t), intent(in)          :: iunits
    type(counter_t), intent(inout)      :: local_counter
    logical, intent(out)                :: boole_t_finished
    real(dp), dimension(3)              :: x_intersection

    if (iper_phi.ne.0) then
        local_counter%phi_0_mappings = local_counter%phi_0_mappings + 1!iper_phi
        if ((local_counter%phi_0_mappings.gt.10) &
            .and.(local_counter%phi_0_mappings.le.3000)) then
            !$omp critical
                write(iunits%pm,*) x
            !$omp end critical
        endif
    endif

    if (x(3).lt.-105d0) then
        boole_t_finished = .true.
       if (local_counter%phi_0_mappings.gt.10) then
        x_intersection = x
        call calc_plane_intersection(x_save,x_intersection,-105d0)
        !$omp critical
            write(iunits%di,*) x, n
        !$omp end critical
       endif
    endif
end subroutine calc_and_write_poincare_mappings_and_divertor_intersections

subroutine calc_plane_intersection(x_save,x,z_plane)

    use constants, only : pi
    
    real(dp), dimension(3), intent(in) :: x_save
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(in) :: z_plane
    real(dp) :: rel_dist_z
    
    rel_dist_z = (z_plane-x_save(3))/(x(3)-x_save(3))
    x(1) = x_save(1) + rel_dist_z*(x(1)-x_save(1))
    if (abs(x(2)-x_save(2)).gt.pi) then
    x(2) = modulo(x_save(2) + 2*pi-abs(x(2)-x_save(2)),2*pi)
    else
    x(2) = x_save(2) + rel_dist_z*(x(2)-x_save(2))
    endif
    x(3) = z_plane
    
    end subroutine calc_plane_intersection

subroutine categorize_lost_particles(ind_tetr,x,pflux,local_counter,t_remain,t_remain_out)

    use tetra_physics_mod, only: tetra_physics
    use boltzmann_types_mod, only: counter_t, poloidal_flux_t

    integer, intent(in) :: ind_tetr
    real(dp), dimension(3), intent(in) :: x
    type(poloidal_flux_t), intent(in) :: pflux
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

subroutine update_local_tetr_moments(moment_specs,local_tetr_moments,ind_tetr,n,start,optional_quantities)

    use boltzmann_types_mod, only: moment_specs_t, start_t
    use gorilla_settings_mod, only: optional_quantities_type

    type(moment_specs_t), intent(in)             :: moment_specs
    type(start_t), intent(in)                    :: start
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

subroutine calc_particle_weights_and_jperp(b,n,z_save,vpar,vperp,ind_tetr,pflux,start)

    use boltzmann_types_mod, only: boltzmann_input_t, poloidal_flux_t, start_t
    use tetra_physics_mod, only: tetra_physics,particle_mass,particle_charge,cm_over_e
    use constants, only: ev2erg
    use volume_integrals_and_sqrt_g_mod, only: sqrt_g
    use supporting_functions_mod, only: bmod_func

    type(boltzmann_input_t), intent(in) :: b
    type(poloidal_flux_t), intent(in) :: pflux
    real(dp), intent(in) :: vpar, vperp
    real(dp), dimension(3), intent(in) :: z_save
    integer, intent(in) :: n,ind_tetr
    type(start_t), intent(inout) :: start
    real(dp) :: local_poloidal_flux, phi_elec_func, temperature
    real(dp) :: r, phi, z

    r = z_save(1)
    phi = z_save(2)
    z = z_save(3)

    if (.not.b%boole_refined_sqrt_g) start%weight = start%weight*r
    if (b%boole_refined_sqrt_g) then
        start%weight(n) = start%weight(n)* (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
                                        &  (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))               
    endif

    if (b%boole_linear_density_simulation.or.b%boole_linear_temperature_simulation) then
        local_poloidal_flux = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*z_save) + &
                                        & cm_over_e*vpar*&
                                        & (tetra_physics(ind_tetr)%h2_1+sum(tetra_physics(ind_tetr)%gh2*z_save))
    endif
    if (b%boole_linear_density_simulation) then
        start%weight(n) = start%weight(n)*(pflux%max*1.1-local_poloidal_flux)/(pflux%max*1.1)
    endif

    if (b%boole_boltzmann_energies) then
        !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts have been added before)
        phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z_save)
        if (.not. b%boole_linear_temperature_simulation) then
            start%weight(n) =start%weight(n)*sqrt(start%energy(n)*ev2erg)/(b%energy_eV*ev2erg)**1.5* &
                        & exp(-(start%energy(n)*ev2erg+particle_charge*phi_elec_func)/(b%energy_eV*ev2erg))
        else
            temperature = b%energy_eV*ev2erg*(pflux%max*1.1-local_poloidal_flux)/(pflux%max*1.1)
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

    if (num_particles.gt.10) then
        if (modulo(kpart,int(num_particles/10)).eq.0) then
            print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
        endif
    else 
        print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
    endif

end subroutine print_progress

subroutine handle_lost_particles(iunit,t_confined, x, n, local_counter, boole_particle_lost)

    use boltzmann_types_mod, only: counter_t

    integer, intent(in) :: iunit, n
    real(dp), intent(in) :: t_confined
    real(dp), dimension(3), intent(in) :: x
    type(counter_t) :: local_counter
    logical :: boole_particle_lost

    !write another if clause (if hole size = minimal .and. particle lost inside .and. boole_cut_out_hole = .true. 
    !(use extra variable m in orbit routine (0 normally, 1 when lost outside, -1 when lost inside))),
    !if m = 1 do as now, if m = -1 select arbitrary newp position and update x, vpar and vperp)
    write(iunit,*) t_confined, x, n
    !$omp critical
    local_counter%lost_particles = 1
    boole_particle_lost = .true.
    !$omp end critical

end subroutine handle_lost_particles

subroutine initialise_loop_variables(b, l, n, v0, start, local_counter,boole_particle_lost,t_step,t_confined, &
    local_tetr_moments,x,v,vpar,vperp)

use boltzmann_types_mod, only: boltzmann_input_t, counter_t, start_t
use constants, only: ev2erg
use tetra_physics_mod, only: particle_mass

type(boltzmann_input_t), intent(in) :: b
type(start_t), intent(in) :: start
integer, intent(in) :: l, n
type(counter_t) :: local_counter
logical :: boole_particle_lost
real(dp) :: t_step, t_confined, pitchpar, v, vpar, vperp, v0
real(dp), dimension(3) :: x
complex(dp), dimension(:,:) :: local_tetr_moments


call set_local_counter_zero(local_counter)
boole_particle_lost = .false.
t_step = b%time_step
t_confined = 0
if (l.eq.1) local_tetr_moments = 0
pitchpar = start%pitch(n)
x = start%x(:,n)
vpar = pitchpar * v0
vperp = sqrt(v0**2-vpar**2)
if (b%boole_boltzmann_energies) then
v = sqrt(start%energy(n)*ev2erg*2/particle_mass)
vpar = pitchpar * v
vperp = sqrt(v**2-vpar**2)
endif

end subroutine initialise_loop_variables

subroutine add_local_tetr_moments_to_output(local_tetr_moments, output, moment_specs)

    use boltzmann_types_mod, only: output_t, moment_specs_t
    use tetra_grid_mod, only: ntetr
    
    complex(dp), dimension(:,:), intent(in) :: local_tetr_moments
    type(moment_specs_t), intent(in) :: moment_specs
    type(output_t) :: output
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

subroutine normalise_prism_moments_and_prism_moments_squared(moment_specs,output,b)

    use boltzmann_types_mod, only: moment_specs_t, output_t, boltzmann_input_t
    
    type(boltzmann_input_t), intent(in) :: b
    type(moment_specs_t), intent(inout) :: moment_specs
    type(output_t), intent(inout) :: output
    integer :: n
    
    
    do n = 1,moment_specs%n_moments
    output%prism_moments(n,:) = output%prism_moments(n,:)/(output%prism_volumes*b%time_step*b%n_particles)
    if (moment_specs%boole_squared_moments) then
    output%prism_moments_squared(n,:) = output%prism_moments_squared(n,:)/ &
                   (output%prism_volumes**2*b%time_step**2*b%n_particles)
    endif
    if (b%boole_refined_sqrt_g) then
    output%prism_moments(n,:) = output%prism_moments(n,:)*output%prism_volumes/output%refined_prism_volumes
    if (moment_specs%boole_squared_moments) then
    output%prism_moments_squared(n,:) = output%prism_moments_squared(n,:)* &
           output%prism_volumes**2/output%refined_prism_volumes**2
    endif
    endif
    enddo
    
end subroutine normalise_prism_moments_and_prism_moments_squared

subroutine set_moment_specifications(moment_specs, boole_squared_moments)

    use gorilla_settings_mod, only: boole_array_optional_quantities
    use tetra_grid_settings_mod, only: grid_size
    use tetra_grid_mod, only: ntetr
    use boltzmann_types_mod, only: moment_specs_t
    
    logical, intent(in) :: boole_squared_moments
    type(moment_specs_t), intent(out) :: moment_specs
    integer :: i, n_prisms
    
    n_prisms = ntetr/3
    
    moment_specs%boole_squared_moments = boole_squared_moments
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

subroutine initialise_output(output, moment_specs)

    use tetra_grid_mod, only: ntetr
    use boltzmann_types_mod, only: moment_specs_t, output_t
    
    type(moment_specs_t), intent(in) :: moment_specs
    type(output_t),intent(out) :: output
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

subroutine fourier_transform_moments(output,moment_specs)

    use constants, only: pi
    use tetra_grid_settings_mod, only: grid_size
    use boltzmann_types_mod, only: moment_specs_t, output_t
    use tetra_grid_mod, only: ntetr
    
    type(output_t), intent(inout)               :: output
    type(moment_specs_t), intent(inout)         :: moment_specs
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

subroutine add_local_counter_to_counter(local_counter,counter)

    use boltzmann_types_mod, only: counter_t
    
    type(counter_t), intent(in) :: local_counter
    type(counter_t), intent(inout) :: counter
    
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

subroutine calc_poloidal_flux(pflux, verts)

    use boltzmann_types_mod, only: poloidal_flux_t
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: ntetr, tetra_grid
    
    real(dp), dimension(:,:), intent(in) :: verts
    integer :: i
    type(poloidal_flux_t), intent(out) :: pflux
    
    pflux%max = 0
    pflux%min = tetra_physics(1)%Aphi1
    do i = 1, ntetr
    pflux%max = max(pflux%max,tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
    & (verts([1,2,3],tetra_grid(i)%ind_knot(4))-verts([1,2,3],tetra_grid(i)%ind_knot(1)))))
    pflux%min = min(pflux%min,tetra_physics(i)%Aphi1)
    enddo
    
end subroutine calc_poloidal_flux

subroutine calc_starting_conditions(b, v0, start, verts)

    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, ntetr
    use find_tetra_mod, only: find_tetra
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_physics_mod, only: coord_system
    use boltzmann_types_mod, only: boltzmann_input_t, start_t
    
    type(boltzmann_input_t), intent(in)                    :: b
    type(start_t)                                          :: start
    real(dp), intent(in)                                   :: v0
    real(dp), dimension(:,:), allocatable, intent(out)     :: verts
    real(dp)                                               :: constant_part_of_weight
    real(dp)                                               :: rand_scalar, vpar, vperp
    real(dp)                                               :: amin, cmin, cmax !amax is set globally
    real(dp), dimension(:), allocatable                    :: rand_vector
    real(dp), dimension(:,:), allocatable                  :: rand_matrix1, rand_matrix2
    integer                                                :: i
    real(dp), dimension(3)                                 :: x
    integer                                                :: ind_tetr_out,iface,seed_inp_unit
    real(dp)                                               :: amax
    integer                                                :: ind_a, ind_b, ind_c
    
    
    !!!!comment out the following section to make starting conditions really random!!!
    
    integer,dimension(:), allocatable                              :: seed
    integer                                                        :: n
    
    print*, 'calc_starting_conditions started'
    open(newunit = seed_inp_unit, file='seed.inp', status='old',action = 'read')
    read(seed_inp_unit,*) n
    allocate(seed(n))
    read(seed_inp_unit,*) seed
    close(seed_inp_unit)
    CALL RANDOM_SEED (PUT=seed)
    deallocate(seed)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate(rand_vector(b%num_particles))
    allocate(rand_matrix1(3,b%num_particles))
    allocate(rand_matrix2(2,b%num_particles))
    
    allocate(start%x(3,b%num_particles))
    allocate(start%pitch(b%num_particles))
    allocate(start%energy(b%num_particles))
    allocate(start%weight(b%num_particles))
    allocate(start%jperp(b%num_particles))
    
    ind_a = 1 !(R in cylindrical coordinates, s in flux coordinates)
    ind_b = 2 !(phi in cylindrical and flux coordinates)
    ind_c = 3 !(z in cylindrical coordinates, theta in flux coordinates)
    
    if (coord_system.eq.2) then
    ind_b = 3
    ind_c = 2
    endif
    
    if (coord_system.eq.1) allocate(verts(size(verts_rphiz(:,1)),size(verts_rphiz(1,:))))
    if (coord_system.eq.2) allocate(verts(size(verts_sthetaphi(:,1)),size(verts_sthetaphi(1,:))))
    if (coord_system.eq.1) verts = verts_rphiz
    if (coord_system.eq.2) verts = verts_sthetaphi
    
    amin = minval(verts(ind_a,:))
    amax = maxval(verts(ind_a,:))
    cmin = minval(verts(ind_c,:))
    cmax = maxval(verts(ind_c,:))
    
    constant_part_of_weight = b%density*(amax-amin)*(cmax-cmin)*2*pi
    
    !compute starting conditions
    if (b%boole_point_source) then
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
    call RANDOM_NUMBER(rand_matrix1)
    start%x(ind_a,:) = amin + (amax - amin)*rand_matrix1(ind_a,:) !r in cylindrical, s in flux coordinates
    start%x(ind_b,:) = 2*pi*rand_matrix1(ind_b,:) !phi in cylindrical and flux coordinates
    start%x(ind_c,:) = cmin + (cmax - cmin)*rand_matrix1(ind_c,:) !z in cylindrical, theta in flux coordinates
    ! start%x(ind_a,:) = (/(214 + i*(216-214)/n_particles, i=1,b%num_particles)/)!r
    ! start%x(ind_b,:) = 0.0d0  !1d-1 !phi in cylindrical and flux coordinates
    ! start%x(ind_c,:) = 12d0 !z in cylindrical, theta in flux coordinates
    endif
    
    call RANDOM_NUMBER(rand_matrix2)
    start%pitch(:) = 2*rand_matrix2(1,:)-1 !pitch parameter
    !start%pitch = 0.7d0 !pitch parameter
    
    if (b%boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
    start%energy = 5*b%energy_eV*rand_matrix2(2,:) !boltzmann energy distribution
    start%weight =  constant_part_of_weight*10/sqrt(pi)*b%energy_eV*ev2erg
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (b%boole_antithetic_variate) then
    start%x(:,1:b%num_particles:2) = start%x(:,2:b%num_particles:2)
    start%pitch(1:b%num_particles:2) = -start%pitch(2:b%num_particles:2)
    start%energy(1:b%num_particles:2) = start%energy(2:b%num_particles:2)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, 'calc_starting_conditions finished'
end subroutine calc_starting_conditions

subroutine calc_collision_coefficients_for_all_tetrahedra(b,c,v0,energy_eV)

    use boltzmann_types_mod, only: collisions_t, boltzmann_input_t
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_physics_mod, only: particle_mass,particle_charge
    use constants, only: echarge,amp
    use tetra_grid_settings_mod, only: grid_size
    use collis_ions, only: collis_init
    
    real(dp), intent(in) :: v0, energy_eV
    type(collisions_t) :: c
    type(boltzmann_input_t) :: b
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
    call collis_init(m0,z0,c%mass_num,c%charge_num,c%dens,c%temp,energy_eV,v0,efcolf,velrat,enrat)
    c%efcolf_mat(:,i) = efcolf
    c%velrat_mat(:,i) = velrat
    c%enrat_mat(:,i) = enrat
    enddo
    
    if (b%boole_precalc_collisions) then
    allocate(c%randcol(b%num_particles,c%randcoli,3))
    call RANDOM_NUMBER(c%randcol)
    !3.464102 = sqrt(12), this creates a random number with zero average and unit variance
    c%randcol(:,:,1:2:3) =  3.464102*(c%randcol(:,:,1:2:3)-.5)
    endif
end subroutine calc_collision_coefficients_for_all_tetrahedra

subroutine carry_out_collisions(b, c, i, n, v0, t, x, vpar, vperp, ind_tetr, iface)

    use boltzmann_types_mod, only: boltzmann_input_t, collisions_t, time_t
    use collis_ions, only: stost
    use find_tetra_mod, only: find_tetra
    
    type(boltzmann_input_t), intent(in) :: b
    type(collisions_t), intent(in) :: c
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
    
    if (b%boole_precalc_collisions) then
    randnum = c%randcol(n,mod(i-1,c%randcoli)+1,:) 
    call stost(efcolf,velrat,enrat,zet,t%step,1,err,(b%time_step-t%confined)*v0,randnum)
    else
    call stost(efcolf,velrat,enrat,zet,t%step,1,err,(b%time_step-t%confined)*v0)
    endif
    
    vpar = zet(5)*zet(4)*v0+vpar_background(1)
    vperp = sqrt(1-zet(5)**2)*zet(4)*v0
    
    !optionally still change particle_mass, particle_charge and cm_over_e, e.g.:
    !particle_charge = particle_charge + echarge
    !particle_mass = particle_mass + ame - amp 
    !cm_over_e = clight*particle_mass/particle_charge
    endif
    
end subroutine carry_out_collisions

end module main_routine_supporting_functions_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module boltzmann_writing_data_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    public :: write_data_to_files, name_files, unlink_files

contains

subroutine write_data_to_files(filenames,output,moment_specs)

    use boltzmann_types_mod, only: filenames_t, boole_writing_data_t, output_t, moment_specs_t

    type(filenames_t), intent(in) :: filenames
    type(output_t), intent(in) :: output
    type(moment_specs_t), intent(in) :: moment_specs
    type(boole_writing_data_t) :: boole_writing_data

    if (boole_writing_data%vertex_indices) &
        call write_vertex_indices(filenames%vertex_indices)

    if (boole_writing_data%vertex_coordinates) &
        call write_vertex_coordinates(filenames%vertex_coordinates)

    if (boole_writing_data%prism_volumes) &
        call write_prism_volumes(filenames%prism_volumes,output%prism_volumes)

    if (boole_writing_data%refined_prism_volumes) &
        call write_refined_prism_volumes(filenames%refined_prism_volumes,output%refined_prism_volumes)

    if (boole_writing_data%boltzmann_density) &
        call write_boltzmann_densities(filenames%boltzmann_density,output%boltzmann_density)

    if (boole_writing_data%electric_potential) &
        call write_electric_potential(filenames%electric_potential,output%electric_potential)

    if (boole_writing_data%moments) then
        if (moment_specs%n_moments.gt.0) then
            call write_moments(filenames%prism_moments,filenames%prism_moments_summed_squares,filenames%tetr_moments,moment_specs, &
                               output%prism_moments,output%prism_moments_squared,output%tetr_moments)
        else
            print*, "Error: moments are not written to file because no moment was computed. Turn computation of moments on in &
                     gorilla.inp."
        endif
    endif

    if (boole_writing_data%fourier_moments) then
        if (moment_specs%n_moments.gt.0) then
            call write_fourier_moments(filenames%fourier_moments,moment_specs, output%moments_in_frequency_space)
        else
            print*, "Error: Fourier moments are not written to file because no moment was computed (and thus also no fourier &
                     moment). Turn computation of moments on in gorilla.inp."
        endif

    endif

end subroutine write_data_to_files

subroutine write_vertex_indices(filename_vertex_indices)

    use tetra_grid_mod, only: ntetr, tetra_grid

    character(len=100) :: filename_vertex_indices
    integer :: vi_unit
    integer :: i

    open(newunit = vi_unit, file = filename_vertex_indices)
    do i=1, ntetr
        write(vi_unit, *) tetra_grid(i)%ind_knot([1, 2, 3, 4])
    end do
    close(vi_unit)

end subroutine write_vertex_indices

subroutine write_vertex_coordinates(filename_vertex_coordinates)

    use tetra_physics_mod, only: coord_system
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, nvert

    character(len=100) :: filename_vertex_coordinates
    integer :: vc_unit
    integer :: i

    open(newunit = vc_unit, file = filename_vertex_coordinates)
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

subroutine write_prism_volumes(filename_prism_volumes,prism_volumes)

    real(dp),dimension(:), intent(in) :: prism_volumes
    character(len=100) :: filename_prism_volumes
    integer :: pv_unit

    open(newunit = pv_unit, file = filename_prism_volumes)
    write(pv_unit,'(ES20.10E4)') prism_volumes
    close(pv_unit)

end subroutine write_prism_volumes

subroutine write_refined_prism_volumes(filename_refined_prism_volumes,refined_prism_volumes)

    real(dp),dimension(:), intent(in) :: refined_prism_volumes
    character(len=100) :: filename_refined_prism_volumes
    integer :: rpv_unit

    open(newunit = rpv_unit, file = filename_refined_prism_volumes)
    write(rpv_unit,'(ES20.10E4)') refined_prism_volumes
    close(rpv_unit)

end subroutine write_refined_prism_volumes

subroutine write_boltzmann_densities(filename_boltzmann_density,boltzmann_density)

    real(dp),dimension(:), intent(in) :: boltzmann_density
    character(len=100) :: filename_boltzmann_density
    integer :: bd_unit

    open(newunit = bd_unit, file = filename_boltzmann_density)
    write(bd_unit,'(ES20.10E4)') boltzmann_density
    close(bd_unit)

end subroutine write_boltzmann_densities

subroutine write_electric_potential(filename_electric_potential,electric_potential)

    real(dp),dimension(:), intent(in) :: electric_potential
    character(len=100) :: filename_electric_potential
    integer :: epv_unit

    open(newunit = epv_unit, file = filename_electric_potential)
    write(epv_unit,'(ES20.10E4)') electric_potential
    close(epv_unit)

end subroutine write_electric_potential

subroutine write_moments(filename_prism_moments, filename_prism_moments_summed_squares, filename_tetr_moments, moment_specs, &
                         prism_moments, prism_moments_squared, tetr_moments)

    use tetra_grid_mod, only: ntetr
    use boltzmann_types_mod, only: moment_specs_t

    complex(dp), dimension(:,:), intent(in) :: tetr_moments, prism_moments, prism_moments_squared
    type(moment_specs_t), intent(in) :: moment_specs
    character(len=100) :: filename_prism_moments, filename_prism_moments_summed_squares, filename_tetr_moments
    integer :: p_moments_unit, pmss_unit, t_moments_unit
    integer :: l, i
    integer :: n_prisms

    open(newunit = p_moments_unit, file = filename_prism_moments)
    open(newunit = pmss_unit, file = filename_prism_moments_summed_squares)
    open(newunit = t_moments_unit, file = filename_tetr_moments)

    n_prisms = ntetr/3

    if (moment_specs%n_moments.gt.0) then
        do l = 1,n_prisms
            do i = 1,moment_specs%n_moments - 1
                write(p_moments_unit,'(2ES20.10E4)',advance="no") real(prism_moments(i,l)), aimag(prism_moments(i,l))
            enddo
                write(p_moments_unit,'(2ES20.10E4)') real(prism_moments(moment_specs%n_moments,l)), &
                                                     aimag(prism_moments(moment_specs%n_moments,l))
        enddo
        if (moment_specs%boole_squared_moments) then
            do l = 1,n_prisms
                do i = 1,moment_specs%n_moments - 1
                    write(pmss_unit,'(2ES20.10E4)',advance="no") real(prism_moments_squared(i,l)), &
                                                                    & aimag(prism_moments_squared(i,l))
                enddo
                    write(pmss_unit,'(2ES20.10E4)') real(prism_moments_squared(moment_specs%n_moments,l)), &
                                                    aimag(prism_moments_squared(moment_specs%n_moments,l))
            enddo
        endif
        do l = 1,ntetr
            do i = 1,moment_specs%n_moments - 1
                write(t_moments_unit,'(2ES20.10E4)',advance="no") real(tetr_moments(i,l)), aimag(tetr_moments(i,l))
            enddo
                write(t_moments_unit,'(2ES20.10E4)') real(tetr_moments(moment_specs%n_moments,l)), &
                                                     aimag(tetr_moments(moment_specs%n_moments,l))
        enddo
    endif

    close(p_moments_unit)
    close(pmss_unit)
    close(t_moments_unit)

end subroutine write_moments

subroutine write_fourier_moments(filename_fourier_moments,moment_specs,moments_in_frequency_space)

    use tetra_grid_settings_mod, only: grid_size
    use boltzmann_types_mod, only: moment_specs_t

    complex(dp), dimension(:,:,:), intent(in) :: moments_in_frequency_space
    type(moment_specs_t), intent(in) :: moment_specs
    character(len=100) :: filename_fourier_moments
    integer :: fm_unit, l, i

    open(newunit = fm_unit, file = filename_fourier_moments)
    do l = 1,moment_specs%n_triangles
        do i = 1,moment_specs%n_fourier_modes-1
            write(fm_unit,'(2ES20.10E4)',advance="no") real(moments_in_frequency_space(1,l,i)), &
                                                       aimag(moments_in_frequency_space(1,l,i))
        enddo
            write(fm_unit,'(2ES20.10E4)') real(moments_in_frequency_space(1,l,moment_specs%n_fourier_modes)), &
                                          aimag(moments_in_frequency_space(1,l,moment_specs%n_fourier_modes))
    enddo
    close(fm_unit)

end subroutine write_fourier_moments

subroutine name_files(filenames)

    use boltzmann_types_mod, only: filenames_t

    type(filenames_t), intent(out) :: filenames

    filenames%exit_times = 'exit_times.dat'
    filenames%remaining_particles = 'remaining_particles.dat'
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

end subroutine name_files

subroutine unlink_files(filenames)

    use boltzmann_types_mod, only: filenames_t

    type(filenames_t) :: filenames

    call unlink(filenames%exit_times)
    call unlink(filenames%remaining_particles)
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

end subroutine unlink_files

subroutine open_files_before_main_loop(iunits)

    use boltzmann_types_mod, only: iunits_t
    
    type(iunits_t) :: iunits
    
    open(newunit = iunits%et, file = 'exit_times.dat')
    open(newunit = iunits%rp, file = 'remaining_particles.dat')
    open(newunit = iunits%pm, file = 'poincare_maps.dat')
    open(newunit = iunits%di, file = 'divertor_intersections.dat')
    
end subroutine open_files_before_main_loop

subroutine close_files(iunits)

    use boltzmann_types_mod, only: iunits_t
    
    type(iunits_t) :: iunits
    
    close(iunits%et)
    close(iunits%rp)
    close(iunits%pm)
    close(iunits%di)
    
end subroutine close_files

end module boltzmann_writing_data_mod