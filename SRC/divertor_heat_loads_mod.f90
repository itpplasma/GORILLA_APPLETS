module divertor_heat_loads_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type iunits_t
    integer :: pm
    integer :: di
    end type iunits_t

    type(iunits_t) :: iunits
   
contains

subroutine read_divertor_heat_loads_inp_into_type

    use boltzmann_types_mod, only: in

    real(dp) :: time_step,energy_eV,n_particles, density, lambda
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_boltzmann_energies, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
               boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
               boole_divertor_intersection, boole_poincare_plot
    integer :: i_integrator_type, seed_option, n_poincare_mappings, n_mappings_ignored

    integer :: dhl_inp_unit
    
    !Namelist for divertor_heat_loads input
    NAMELIST /divertor_heat_loads_nml/ time_step,energy_eV,n_particles, lambda, boole_divertor_intersection, boole_poincare_plot, &
    & n_poincare_mappings,n_mappings_ignored,boole_squared_moments,boole_point_source,boole_collisions, boole_precalc_collisions, &
    & density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, boole_antithetic_variate, &
    & boole_linear_temperature_simulation,i_integrator_type,seed_option, boole_write_vertex_indices, &
    & boole_write_vertex_coordinates, boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
    & boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data

    open(newunit = dhl_inp_unit, file='divertor_heat_loads.inp', status='unknown')
    read(dhl_inp_unit,nml=divertor_heat_loads_nml)
    close(dhl_inp_unit)

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
    in% boole_divertor_intersection =  boole_divertor_intersection
    in%boole_poincare_plot = boole_poincare_plot
    in%n_poincare_mappings = n_poincare_mappings
    in%n_mappings_ignored = n_mappings_ignored
    in%lambda = lambda

    print *,'GORILLA: Loaded input data from boltzmann.inp'

end subroutine read_divertor_heat_loads_inp_into_type

subroutine calc_divertor_heat_loads

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg,pi,echarge
    use tetra_physics_mod, only: particle_mass,particle_charge
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_set_num_threads
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_grid_mod, only: ntetr, verts_rphiz
    use gorilla_settings_mod, only: ispecies
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use boltzmann_types_mod, only: moment_specs, counter_t, counter, c, in, time_t, particle_status_t
    use write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use calc_boltzmann_supporting_functions_mod, only: print_progress, handle_lost_particles, initialise_loop_variables, &
    add_local_tetr_moments_to_output, normalise_prism_moments_and_prism_moments_squared, set_moment_specifications, &
    initialise_output, fourier_transform_moments, add_local_counter_to_counter, get_ipert, calc_poloidal_flux, &
    calc_collision_coefficients_for_all_tetrahedra, carry_out_collisions, initialize_exit_data, update_exit_data, &
    set_seed_for_random_numbers

    integer :: kpart,i,n,l,p,ind_tetr,iface,iantithetic
    real(dp) :: v0,vpar,vperp
    real(dp), dimension(3) :: x
    complex(dp), dimension(:,:), allocatable :: local_tetr_moments
    type(counter_t) :: local_counter
    type(time_t) :: t
    type(particle_status_t) :: particle_status

    call set_seed_for_random_numbers
    call read_divertor_heat_loads_inp_into_type
    call get_ipert()
    call initialize_gorilla(i_option,ipert)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete this again afterwards !!!!!!!!!!!!!!!!!!!!!!!
    if (ispecies.eq.4) particle_charge = 15*echarge
    print*, 'particle charge number = ', particle_charge/echarge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    v0=sqrt(2.0_dp*in%energy_eV*ev2erg/particle_mass)
    kpart = 0
    iantithetic = 1
    if (in%boole_antithetic_variate) iantithetic = 2

    call set_moment_specifications
    call initialise_output
    call calc_square_root_g
    call calc_volume_integrals(in%boole_boltzmann_energies,in%boole_refined_sqrt_g, in%density, in%energy_eV)
    call calc_starting_conditions
    call initialize_exit_data
    call calc_poloidal_flux(verts_rphiz)
    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra(v0)
    call give_file_names
    call unlink_files
    call open_files

    allocate(local_tetr_moments(moment_specs%n_moments,ntetr))
    if (in%i_integrator_type.eq.2) print*, 'Error: i_integratpr_type set to 2, this module only works with &
                                    & i_integrator_type set to 1'

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(counter, kpart,v0, in, c, iantithetic) &
    !$OMP& PRIVATE(p,l,n,i,x,vpar,vperp,t,ind_tetr,iface,local_tetr_moments,local_counter,particle_status)
    print*, 'get number of threads', omp_get_num_threads()
    !$OMP DO

    !Loop over particles
    do p = 1,in%num_particles/iantithetic
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
                ! if (i.eq.1) then
                !     particle_status%initialized = .false.
                ! endif

                if (in%boole_collisions) then
                    call carry_out_collisions(i, n, v0, t, x, vpar,vperp,ind_tetr, iface)
                    t%step = t%step/v0 !in carry_out_collisions, t%step is initiated as a length, so you need to divide by v0
                endif

                call orbit_timestep_dhl(x,vpar,vperp,t%step,particle_status,ind_tetr, &
                                        iface,n,local_tetr_moments,local_counter,t%remain)

                t%confined = t%confined + t%step - t%remain

                if (ind_tetr.eq.-1) then
                    call handle_lost_particles(local_counter, particle_status%lost)
                    exit
                endif
                if (particle_status%exit) exit
            enddo

            !$omp critical
            counter%integration_steps = counter%integration_steps + i
            c%maxcol = max(dble(i)/dble(c%randcoli),c%maxcol)
            call add_local_counter_to_counter(local_counter)
            !$omp end critical
            call update_exit_data(particle_status%lost,t%confined,x,vpar,vperp,i,n,local_counter%phi_0_mappings)
        enddo
        !$omp critical
        call add_local_tetr_moments_to_output(local_tetr_moments)
        !$omp end critical
    enddo !n
    !$OMP END DO
    !$OMP END PARALLEL

    call normalise_prism_moments_and_prism_moments_squared
    if (moment_specs%n_moments.gt.0) call fourier_transform_moments
    call close_files
    call write_data_to_files
    !call create_magnetic_field_file

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
         print*, 'number of particles left through the outside = ', counter%lost_outside
         print*, 'number of particles left through the inside = ', counter%lost_inside
    endif

end subroutine calc_divertor_heat_loads

subroutine orbit_timestep_dhl(x,vpar,vperp,t_step,particle_status,ind_tetr,iface, n,local_tetr_moments, local_counter, t_remain_out)

    use pusher_tetra_rk_mod, only: pusher_tetra_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly
    use tetra_physics_mod, only: tetra_physics
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: vperp_func, bmod_func
    use find_tetra_mod, only: find_tetra
    use boltzmann_types_mod, only: counter_t, particle_status_t, in
    use tetra_grid_settings_mod, only: grid_kind
    use orbit_timestep_gorilla_supporting_functions_mod, only: categorize_lost_particles, update_local_tetr_moments, &
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
    logical                                      :: boole_t_finished
    integer                                      :: ind_tetr_save,iper_phi,n
    type(optional_quantities_type)               :: optional_quantities
    
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
    ! if (n.eq.1) then
    !     print*, 'hello', x, bmod_func(z_save,ind_tetr)
    !     stop
    ! endif
    !open(581,file = 'banana_orbit.dat')

    do !Loop for tetrahedron pushings until t_step is reached
        local_counter%tetr_pushings = local_counter%tetr_pushings +1

        if(ind_tetr.eq.-1) then
            if((grid_kind.eq.2).or.(grid_kind.eq.3)) then
                call categorize_lost_particles(ind_tetr_save,z_save,local_counter,t_remain,t_remain_out)
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
        !write(581,*) x

        t_remain = t_remain - t_pass

        !call calc_and_write_poincare_mappings_and_divertor_intersections(x_save,x,n,iper_phi,local_counter,particle_status)
        if(in%boole_poincare_plot)          call calc_and_write_poincare_mappings(x,iper_phi,local_counter,particle_status)
        if(in%boole_divertor_intersection)  call calc_and_write_divertor_intersections(x_save,x,n,local_counter,particle_status)
        call update_local_tetr_moments(local_tetr_moments,ind_tetr_save,n,optional_quantities)
        if(particle_status%exit.or.boole_t_finished) then
            if( present(t_remain_out)) t_remain_out = t_remain
            exit
        endif

    enddo !Loop for tetrahedron pushings
    !close(581)

    vperp = vperp_func(z_save,perpinv,ind_tetr_save) !Compute vperp from position

end subroutine orbit_timestep_dhl

subroutine calc_and_write_poincare_mappings(x,iper_phi,local_counter,particle_status)

    use boltzmann_types_mod, only: counter_t, particle_status_t, in

    real(dp), dimension(3), intent(inout)  :: x
    integer, intent(in)                    :: iper_phi
    type(counter_t), intent(inout)         :: local_counter
    type(particle_status_t), intent(inout) :: particle_status

    if (iper_phi.ne.0) then
        local_counter%phi_0_mappings = local_counter%phi_0_mappings + 1!iper_phi
        if ((local_counter%phi_0_mappings.ge.in%n_mappings_ignored) &
            .and.(local_counter%phi_0_mappings.le.in%n_poincare_mappings)) then
            !$omp critical
                write(iunits%pm,*) x
            !$omp end critical
        endif
        if (local_counter%phi_0_mappings.ge.in%n_poincare_mappings) particle_status%exit = .true.
    endif

end subroutine calc_and_write_poincare_mappings

subroutine calc_and_write_divertor_intersections(x_save,x,n,local_counter,particle_status)

    use boltzmann_types_mod, only: counter_t, particle_status_t, in

    real(dp), dimension(3), intent(in)     :: x_save
    real(dp), dimension(3), intent(inout)  :: x
    integer, intent(in)                    :: n
    type(counter_t), intent(inout)         :: local_counter
    type(particle_status_t), intent(inout) :: particle_status
    real(dp), dimension(3)                 :: x_intersection
    real(dp)                               :: z_divertor_plate = -105.0_dp

    if (x(3).lt.z_divertor_plate) then
        particle_status%exit = .true.
       if (local_counter%phi_0_mappings.gt.in%n_mappings_ignored) then
            x_intersection = x
            call calc_plane_intersection(x_save,x_intersection,z_divertor_plate)
            x = x_intersection
            !$omp critical
                write(iunits%di,*) x, n
            !$omp end critical
       endif
    endif

end subroutine calc_and_write_divertor_intersections

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


subroutine open_files
    
    open(newunit = iunits%pm, file = 'poincare_maps.dat')
    open(newunit = iunits%di, file = 'divertor_intersections.dat')
    
end subroutine open_files

subroutine close_files
    
    close(iunits%pm)
    close(iunits%di)
    
end subroutine close_files

subroutine calc_starting_conditions

    use boltzmann_types_mod, only: in
    
    real(dp), dimension(:,:), allocatable  :: rand_matrix

    call set_coordinate_limits

    allocate(rand_matrix(5,in%num_particles))
    call RANDOM_NUMBER(rand_matrix)

    call allocate_start_type
    call set_start_type(rand_matrix)

end subroutine calc_starting_conditions

subroutine set_start_type(rand_matrix)

    use boltzmann_types_mod, only: in, start, g
    use constants, only: pi, ev2erg

    real(dp), dimension(:,:), intent(in) :: rand_matrix
    integer                              :: i
    real(dp)                             :: constant

    if (in%num_particles.gt.1) start%x(1,:) = (/(214 + (i-1)*(216-214)/(in%n_particles-1), i=1,in%num_particles)/)!r
    if (in%num_particles.eq.1) start%x(1,:) = 214
    start%x(2,:) = 0.0_dp  !phi
    start%x(3,:) = 12.0_dp !z
    !start%pitch(:) = 2*rand_matrix(4,:)-1

    constant = (1-in%lambda**2)*start%x(1,1)
    do i = 1,in%num_particles
        start%pitch(i) = sqrt(1-constant/start%x(1,i))
    enddo

    start%weight = in%density*(g%amax-g%amin)*(g%cmax-g%cmin)*2*pi
    start%energy = in%energy_eV

    if (in%boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start%weight =  start%weight*10/sqrt(pi)*in%energy_eV*ev2erg
        start%energy = 5*in%energy_eV*rand_matrix(5,:)
    endif
    
    if (in%boole_antithetic_variate) then
        start%x(:,1:in%num_particles:2) = start%x(:,2:in%num_particles:2)
        start%pitch(1:in%num_particles:2) = -start%pitch(2:in%num_particles:2)
        start%energy(1:in%num_particles:2) = start%energy(2:in%num_particles:2)
    endif

end subroutine set_start_type

subroutine set_coordinate_limits

    use tetra_grid_mod, only: verts_rphiz
    use boltzmann_types_mod, only: g

    g%amin = minval(verts_rphiz(1,:)) !r coordinate
    g%amax = maxval(verts_rphiz(1,:))
    g%cmin = minval(verts_rphiz(3,:)) !z coordinate
    g%cmax = maxval(verts_rphiz(3,:))

end subroutine set_coordinate_limits

subroutine allocate_start_type

    use boltzmann_types_mod, only: start, in

    allocate(start%x(3,in%num_particles))
    allocate(start%pitch(in%num_particles))
    allocate(start%energy(in%num_particles))
    allocate(start%weight(in%num_particles))
    allocate(start%jperp(in%num_particles))

end subroutine allocate_start_type

subroutine create_magnetic_field_file
    use tetra_grid_settings_mod, only: n1,n2,n3
    use boltzmann_types_mod, only: g
    use constants, only: pi
    use tetra_physics_mod, only: vector_potential_rphiz
    use field_mod, only: ipert

    character(len=100) :: filename
    real(dp) :: A(3), B(3), bmod
    real(dp), dimension(3,2) :: limits !min/max for the three coordinates
    integer, dimension(3) :: n_points !number of vetices per direction
    logical :: is_periodic(3) !describes periodicity in three dimensions
    integer :: unit, i, j, k, n
    real(dp), dimension(:), allocatable :: R, phi, Z

    filename = "field_from_gorilla"

    limits(1,:) = (/g%amin,g%amax/) !R
    limits(2,:) = (/0.0_dp,2*pi/) !phi
    limits(3,:) = (/g%cmin,g%amax/) !Z
    n_points = (/n1, n2, n3/)
    is_periodic = (/.false.,.true.,.false./)

    allocate(R(n1), phi(n2), Z(n3))
    R = linspace(limits(1,1), limits(1,2), n1)
    phi = linspace(limits(2,1), limits(2,2), n2)
    Z = linspace(limits(3,1), limits(3,2), n3)

    open(newunit = unit, file = filename)
    write(unit,*) n_points
    write(unit,*) limits
    write(unit,*) is_periodic

    do i = 1, n1
        do j = 1, n2
            do k = 1, n3
                call vector_potential_rphiz(R(i), phi(j), Z(k),ipert,1.0_dp,A(1),A(2),A(3),B(1),B(2),B(3),bmod)
                write(unit,*) A, B
            end do
        end do
    end do

    close(unit)

end subroutine create_magnetic_field_file

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

end module divertor_heat_loads_mod