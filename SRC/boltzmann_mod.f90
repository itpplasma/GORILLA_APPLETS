module boltzmann_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    public :: calc_boltzmann

    private

    real(dp), dimension(:,:), allocatable :: verts
    real(dp), dimension(:), allocatable :: weights, J_perp, temperature_vector
    integer :: et_unit, rp_unit, pm_unit, di_unit
   
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
    use tetra_grid_settings_mod, only: n_field_periods, grid_size
    use tetra_grid_mod, only: ntetr, nvert, verts_rphiz, tetra_grid
    use gorilla_settings_mod, only: ispecies
    use collis_ions, only: stost
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use boltzmann_types_mod, only: filenames_t, output_t, moment_specs_t, counter_t, poloidal_flux_t, collisions_t,boltzmann_input_t
    use boltzmann_writing_data_mod, only: write_data_to_files, name_files, unlink_files

    real(dp), dimension(:,:), allocatable :: start_pos_pitch_mat
    real(dp) :: v0,pitchpar,vpar,vperp,t_remain,t_confined, v, maxcol
    integer :: kpart,i,j,n,l,m,k,p,ind_tetr,iface,ierr,err,iantithetic, num_background_species, inorout
    integer :: i_part
    real(dp), dimension(3) :: x_rand_beg,x,randnum
    logical :: boole_initialized,boole_particle_lost
    real(dp) :: dtau, dphi,dtaumin, t_step
    real(dp), dimension(5) :: z, zet
    Character(LEN=50) :: format_moments, format_fourier_moments
    complex(dp), dimension(:,:), allocatable :: local_tetr_moments
    real(dp) :: m0,z0,hamiltonian_time
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat,vpar_background
    integer :: Te_unit, Ti_unit, ne_unit, t_moments_unit, ipert_unit
    type(filenames_t) :: filenames
    type(output_t) :: output
    type(moment_specs_t) :: moment_specs
    type(counter_t) :: counter, local_counter
    real(dp), dimension(:,:), allocatable :: sqrt_g
    type(poloidal_flux_t) :: poloidal_flux
    type(collisions_t) :: c
    integer :: n_prisms
    real(dp), dimension(:,:,:), allocatable :: randcol
    integer :: randcoli = int(1.0d5)
    type(boltzmann_input_t) :: b

    call read_boltzmann_inp_into_type(b)
    call get_ipert()
    call initialize_gorilla(i_option,ipert)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete this again afterwards !!!!!!!!!!!!!!!!!!!!!!!
    if (ispecies.eq.4) particle_charge = 15*echarge
    print*, 'particle charge number = ', particle_charge/echarge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    v0=sqrt(2.d0*b%energy_eV*ev2erg/particle_mass)
    n_prisms = ntetr/3
    allocate(weights(b%num_particles))
    allocate(J_perp(b%num_particles))
    allocate(temperature_vector(b%num_particles))
    if (b%boole_precalc_collisions) allocate(randcol(b%num_particles,randcoli,3))
    temperature_vector = 0
    kpart = 0
    maxcol = 0
    iantithetic = 1
    if (b%boole_antithetic_variate) iantithetic = 2

    call set_moment_specifications(moment_specs, b%boole_squared_moments)
    call initialise_output(output, moment_specs)
    call calc_square_root_g(sqrt_g)
    call calc_volume_integrals(b%boole_boltzmann_energies,b%boole_refined_sqrt_g, b%density, b%energy_eV,output)
    call calc_starting_conditions(b,v0,start_pos_pitch_mat,randcol)
    call initialize_poloidal_flux(poloidal_flux,b%num_particles)
    call prepare_collisions(c,v0,b%energy_eV)
    call name_files(filenames)
    call unlink_files(filenames)
    call open_files_before_main_loop

    allocate(local_tetr_moments(moment_specs%n_moments,ntetr))
    if (b%boole_collisions) allocate(efcolf(num_background_species),velrat(num_background_species),enrat(num_background_species))

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(kpart,v0,sqrt_g, output, &
    !$OMP& dtau,dtaumin,counter, &
    !$OMP& start_pos_pitch_mat,et_unit, moment_specs, poloidal_flux, b, c,&
    !$OMP& tetra_grid,iantithetic,tetra_physics, di_unit, pm_unit, rp_unit, &
    !$OMP& num_background_species,randcol,randcoli,maxcol) &
    !$OMP& FIRSTPRIVATE(particle_mass, particle_charge) &
    !$OMP& PRIVATE(p,l,n,boole_particle_lost,x_rand_beg,x,pitchpar,vpar,vperp,boole_initialized,t_step,err,zet, &
    !$OMP& ind_tetr,iface,t_remain,t_confined,z,ierr, v, local_tetr_moments,hamiltonian_time, &
    !$OMP& m0,z0,i,efcolf,velrat,enrat,vpar_background,inorout,randnum,local_counter)
    print*, 'get number of threads', omp_get_num_threads()
    !$OMP DO

    !Loop over particles
    do p = 1,b%num_particles/iantithetic
        do l = 1,iantithetic

            n = (p-1)*iantithetic+l
            !$omp critical
            !Counter for particles
            kpart = kpart+1 !in general not equal to n becuase of parallelisation
            boole_particle_lost = .false.
if (b%num_particles.gt.10) then
if (modulo(kpart,int(b%num_particles/10)).eq.0) then
    print *, kpart, ' / ', b%num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
endif
else 
print *, kpart, ' / ', b%num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
endif
            !$omp end critical

            call set_local_counter_zero(local_counter) 

            t_step = b%time_step
            t_confined = 0
            if (l.eq.1) local_tetr_moments = 0

            !You need x_rand_beg(1,3), pitchpar(1) (between -1 and 1), energy is already given
            x_rand_beg = start_pos_pitch_mat(1:3,n)
            pitchpar = start_pos_pitch_mat(4,n)
            ! print*, start_pos_pitch_mat(5,n)
            ! if (n.eq.972) print*, x_rand_beg,start_pos_pitch_mat(5,n),pitchpar

            if (b%i_integrator_type.eq.2) print*, 'Error: i_integratpr_type set to 2, this module only works with &
                                                & i_integrator_type set to 1'
            x = x_rand_beg
            vpar = pitchpar * v0
            vperp = sqrt(v0**2-vpar**2)
            if (b%boole_boltzmann_energies) then
                v = sqrt(start_pos_pitch_mat(5,n)*ev2erg*2/particle_mass)
                vpar = pitchpar * v
                vperp = sqrt(v**2-vpar**2)
            endif

            i = 0
            do while (t_confined.lt.b%time_step)
                i = i+1
                !Orbit integration
                if (i.eq.1) then
                    boole_initialized = .false.
                endif
                if (b%boole_collisions) then
                    if (i.eq.1) call find_tetra(x,vpar,vperp,ind_tetr,iface)
                    if (.not.(ind_tetr.eq.-1)) then
                        efcolf = c%efcolf_mat(:,ind_tetr)
                        velrat = c%velrat_mat(:,ind_tetr)
                        enrat = c%enrat_mat(:,ind_tetr)
                        vpar_background = c%vpar_mat(:,ind_tetr)
                        !print*, vpar_background

                        vpar = vpar - vpar_background(1)
                        !since vpar_background actually has num_background_particles entries, consider giving it as an extra
                        !optional input variable to stost, before randnum (maybe also check if radnum could then be set by 
                        !randnum = variable eve if vpar_background is not set and other variables are not set by name indexing)
                        !since it came up when writing these lines: replace expressions like
                        !"verts(size(verts_rphiz(:,1)),size(verts_rphiz(1,:)))" with "3,nvert"
                        zet(1:3) = x !spatial position
                        zet(4) = sqrt(vpar**2+vperp**2)/v0 !normalized velocity module 
                        zet(5) = vpar/sqrt(vpar**2+vperp**2) !pitch parameter

                        if (b%boole_precalc_collisions) then
                            randnum = randcol(n,mod(i-1,randcoli)+1,:) 
                            call stost(efcolf,velrat,enrat,zet,t_step,1,err,(b%time_step-t_confined)*v0,randnum)
                        else
                            call stost(efcolf,velrat,enrat,zet,t_step,1,err,(b%time_step-t_confined)*v0)
                        endif

                        t_step = t_step/v0
                        x = zet(1:3)
                        vpar = zet(5)*zet(4)*v0+vpar_background(1)
                        vperp = sqrt(1-zet(5)**2)*zet(4)*v0

                        !optionally still change particle_mass, particle_charge and cm_over_e, e.g.:
                        !particle_charge = particle_charge + echarge
                        !particle_mass = particle_mass + ame - amp 
                        !cm_over_e = clight*particle_mass/particle_charge
                    endif
                endif

                call orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface,n,v, sqrt_g,&
                            & start_pos_pitch_mat,local_tetr_moments, moment_specs, b, poloidal_flux, hamiltonian_time,inorout, &
                              local_counter,t_remain)

                t_confined = t_confined + t_step - t_remain
                !Lost particle handling
                if(ind_tetr.eq.-1) then
!write another if clause (if hole size = minimal .and. particle lost inside .and. boole_cut_out_hole = .true. 
!(use extra variable m in orbit routine (0 normally, 1 when lost outside, -1 when lost inside))),
!if m = 1 do as now, if m = -1 select arbitrary newp position and update x, vpar and vperp)
                    write(et_unit,*) t_confined, x, n
                    !print*, t_confined, x
                    !$omp critical
                    local_counter%lost_particles = 1
                    boole_particle_lost = .true.
                    !$omp end critical
                    exit
                endif

                v = sqrt(vpar**2+vperp**2)
            enddo
            !$omp critical
            counter%integration_steps = counter%integration_steps + i
            maxcol = max(dble(i)/dble(randcoli),maxcol)
            call add_local_counter_to_counter(local_counter,counter)
            !$omp end critical
            if (t_confined.eq.b%time_step) then
                write(rp_unit,*) x, v, vpar, vperp, i, n
            endif
        enddo

        !$omp critical
        call add_local_tetr_moments_to_output(local_tetr_moments, output, moment_specs)
        !$omp end critical
    enddo !n
    !$OMP END DO
    !$OMP END PARALLEL

    call normalise_prism_moments_and_prism_moments_squared(moment_specs,output,b)
    if (moment_specs%n_moments.gt.0) call fourier_transform_moments(output, moment_specs)
    call close_files
    call write_data_to_files(filenames,output,moment_specs)

    if (b%boole_precalc_collisions) print*, "maxcol = ", maxcol
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
    print*, 'number of particles left through the outside = ', counter%lost_outside
    print*, 'number of particles left through the inside = ', counter%lost_inside

end subroutine calc_boltzmann

subroutine orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, n,v,sqrt_g,start_pos_pitch_mat, &
                                          & local_tetr_moments,moment_specs, b, poloidal_flux, hamiltonian_time, inorout, &
                                            local_counter, t_remain_out)

    use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
    use tetra_physics_poly_precomp_mod , only: make_precomp_poly_perpinv, initialize_boole_precomp_poly_perpinv, &
        & alloc_precomp_poly_perpinv
    use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass,cm_over_e
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type, boole_array_optional_quantities
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: bmod_func, vperp_func
    use find_tetra_mod, only: find_tetra
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: tetra_grid, ntetr
    use boltzmann_types_mod, only: moment_specs_t, counter_t, poloidal_flux_t, boltzmann_input_t

    type(moment_specs_t), intent(in)        :: moment_specs
    type(counter_t), intent(inout)          :: local_counter
    type(poloidal_flux_t)                   :: poloidal_flux
    type(boltzmann_input_t)                 :: b
    real(dp), dimension(3), intent(inout)   :: x
    real(dp), intent(inout)                 :: vpar,vperp
    real(dp), intent(in)                    :: t_step
    logical, intent(inout)                  :: boole_initialized
    integer, intent(inout)                  :: ind_tetr,iface
    real(dp), intent(out), optional         :: t_remain_out
    real(dp), dimension(3)                  :: z_save, x_save
    real(dp)                                :: vperp2,t_remain,t_pass,vpar_save, v, hamiltonian_time, aphi
    logical                                 :: boole_t_finished
    integer                                 :: ind_tetr_save,iper_phi,k, m, n, inorout
    real(dp)                                :: perpinv,perpinv2, speed, r, z, phi, B_field, phi_elec_func
    real(dp), dimension(:,:)                :: start_pos_pitch_mat
    type(optional_quantities_type)          :: optional_quantities
    complex(dp), dimension(:,:)             :: local_tetr_moments
    integer                                 :: single_particle_counter_phi0_mappings
    real(dp), dimension(:,:)                :: sqrt_g
    

    !If orbit_timestep is called for the first time without grid position
    if(.not.boole_initialized) then

        !Check coordinate domain (optionally perform modulo operation)
        call check_coordinate_domain(x)

        !Find tetrahedron index and face index for position x
        call find_tetra(x,vpar,vperp,ind_tetr,iface)

        !If particle doesn't lie inside any tetrahedron
        if(ind_tetr.eq.-1) then
            t_remain_out = t_step
            return
        endif

        r = x(1) - verts(1, tetra_grid(ind_tetr)%ind_knot(1))
        z = x(3) - verts(3,tetra_grid(ind_tetr)%ind_knot(1))
        phi = x(2) - verts(2,tetra_grid(ind_tetr)%ind_knot(1))

        if (b%boole_refined_sqrt_g) then
            weights(n) = weights(n)* &
                                    &  (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
                                    &  (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))
        endif

        if (b%boole_linear_density_simulation.or.b%boole_linear_temperature_simulation) then
            poloidal_flux%particle(n) = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*(/r,phi,z/)) + &
                                             & cm_over_e*vpar*&
                                             & (tetra_physics(ind_tetr)%h2_1+sum(tetra_physics(ind_tetr)%gh2*(/r,phi,z/)))
        endif
        if (b%boole_linear_density_simulation) then
            weights(n) = weights(n)*(poloidal_flux%max*1.1-poloidal_flux%particle(n))/(poloidal_flux%max*1.1)
        endif
      
        if (b%boole_boltzmann_energies) then
            !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts have been added before)
            phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*(/r,phi,z/))
            if (.not. b%boole_linear_temperature_simulation) then
                weights(n) = weights(n)*sqrt(start_pos_pitch_mat(5,n)*ev2erg)/(b%energy_eV*ev2erg)**1.5* &
                            & exp(-(start_pos_pitch_mat(5,n)*ev2erg+particle_charge*phi_elec_func)/(b%energy_eV*ev2erg))
            else
                temperature_vector(n) = b%energy_eV*ev2erg*(poloidal_flux%max*1.1-poloidal_flux%particle(n))/(poloidal_flux%max*1.1)
                weights(n) = weights(n)*sqrt(start_pos_pitch_mat(5,n)*ev2erg)/temperature_vector(n)**1.5* &
                & exp(-(start_pos_pitch_mat(5,n)*ev2erg+particle_charge*phi_elec_func)/temperature_vector(n))
            endif
        endif

        !compute J_perp for perpendicular pressure
        B_field = tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*(/r,phi,z/))
        J_perp(n) = particle_mass*v**2*(1-start_pos_pitch_mat(4,n)**2)*cm_over_e/(2*B_field)*(-1)
        !-1 because of negative gyrophase

        boole_initialized = .true.
    endif
          
    !Exit the subroutine after initialization, if time step equals zero
    if(t_step.eq.0.d0) return

    inorout = 0

    !Squared perpendicular velocity
    vperp2 = vperp**2

    !Compute relative particle position
    z_save = x-tetra_physics(ind_tetr)%x1

    !Compute perpendicular invariant of particle
    perpinv=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)
    perpinv2 = perpinv**2
             
    !Initialize constants of motion in particle-private module
    select case(ipusher)
        case(1)
            call initialize_const_motion_rk(perpinv,perpinv2)
        case(2)
            call initialize_const_motion_poly(perpinv,perpinv2)
    end select        
!
!             !NOT FULLY IMPLEMENTED YET: Precompute quantities dependent on perpinv
!             call alloc_precomp_poly_perpinv(1,ntetr)
!             call initialize_boole_precomp_poly_perpinv()
!             call make_precomp_poly_perpinv(perpinv,perpinv2)
!
    !Integrate particle orbit for given time step
    t_remain = t_step

    !Logical for handling time integration
    boole_t_finished = .false.

    single_particle_counter_phi0_mappings = 0

    local_counter%tetr_pushings = local_counter%tetr_pushings -1 !set tetr_pushings to -1 because when entering the loop it will 
    !go back to one without pushing

    !Loop for tetrahedron pushings until t_step is reached
    do
        local_counter%tetr_pushings = local_counter%tetr_pushings +1
        !Domain Boundary
        if(ind_tetr.eq.-1) then
            !print *, 'WARNING: Particle lost.'
            aphi = tetra_physics(ind_tetr_save)%Aphi1+sum(tetra_physics(ind_tetr_save)%gAphi*z_save)
            !$omp critical
            if (abs(aphi-poloidal_flux%max).lt.abs(aphi-poloidal_flux%min)) then
                local_counter%lost_outside = local_counter%lost_outside + 1
                inorout = 1
            endif
            if (abs(aphi-poloidal_flux%max).gt.abs(aphi-poloidal_flux%min)) then
                local_counter%lost_inside = local_counter%lost_inside + 1
                inorout = -1
            endif
            !$omp end critical
            if( present(t_remain_out)) then
                t_remain_out = t_remain
            endif
            exit
        endif
          
        !Save the tetrahedron index for computation of vperp in the last step
        ind_tetr_save = ind_tetr

        !Save vpar for the computation of the parallel adiabatic invariant
        vpar_save = vpar

        !t_remain (in) ... remaining time until t_step is finished
        !t_pass (out) ... time to pass the tetrahdron

        !Save x for the computation of the current
        x_save = x

        !Calculate trajectory
        select case(ipusher)
            case(1)
                call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper_phi)
            case(2)
                call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_remain,&
                                                    & t_pass,boole_t_finished,iper_phi,optional_quantities)
        end select

        t_remain = t_remain - t_pass
        if (iper_phi.ne.0) then
            local_counter%phi_0_mappings = local_counter%phi_0_mappings + 1!iper_phi
            single_particle_counter_phi0_mappings = single_particle_counter_phi0_mappings + 1
            if ((local_counter%phi_0_mappings.gt.10) &
                .and.(local_counter%phi_0_mappings.le.3000)) then
                !$omp critical
                    write(pm_unit,*) x
                !$omp end critical
            endif
        endif

        if (x(3).lt.-105d0) then
            boole_t_finished = .true.
           if (single_particle_counter_phi0_mappings.gt.10) then
            call calc_plane_intersection(x_save,x,-105d0)
            !$omp critical
                write(di_unit,*) x, n
            !$omp end critical
           endif
        endif

        hamiltonian_time = 0


        do m = 1,moment_specs%n_moments
            select case(moment_specs%moments_selector(m))
                case(1)
                    local_tetr_moments(m,ind_tetr_save) = local_tetr_moments(m,ind_tetr_save) + &
                                                                    & weights(n)*optional_quantities%t_hamiltonian!* &
                                                                    !& (exp(2*(0,1)*x(2))+exp(3*(0,1)*x(2)))
                    hamiltonian_time = optional_quantities%t_hamiltonian
                case(2)
                    local_tetr_moments(m,ind_tetr_save) = local_tetr_moments(m,ind_tetr_save) + &
                                                                    & weights(n)*optional_quantities%gyrophase*J_perp(n)
                case(3)
                    local_tetr_moments(m,ind_tetr_save) = local_tetr_moments(m,ind_tetr_save) + &
                                                                    & weights(n)*optional_quantities%vpar_int
                case(4)
                    local_tetr_moments(m,ind_tetr_save) = local_tetr_moments(m,ind_tetr_save) + &
                                                                    & weights(n)*optional_quantities%vpar2_int
            end select
        enddo

        !Orbit stops within cell, because "flight"-time t_step has finished
        if(boole_t_finished) then
            if( present(t_remain_out)) then
                t_remain_out = t_remain
            endif
            exit
        endif

    enddo !Loop for tetrahedron pushings

    !Compute vperp from position
    vperp = vperp_func(z_save,perpinv,ind_tetr_save)
        
!             !NOT FULLY IMPLEMENTED YET: Deallocate precomputed quantities dependent on perpinv
!             call alloc_precomp_poly_perpinv(2,ntetr)

end subroutine orbit_timestep_gorilla_boltzmann

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

subroutine open_files_before_main_loop

    open(newunit = et_unit, file = 'exit_times.dat')
    open(newunit = rp_unit, file = 'remaining_particles.dat')
    open(newunit = pm_unit, file = 'poincare_maps.dat')
    open(newunit = di_unit, file = 'divertor_intersections.dat')

end subroutine open_files_before_main_loop

subroutine close_files

    close(et_unit)
    close(rp_unit)
    close(pm_unit)
    close(di_unit)

end subroutine close_files

subroutine set_local_counter_zero(counter)

    use boltzmann_types_mod, only: counter_t

    type(counter_t), intent(inout) :: counter

    counter%lost_particles = 0
    counter%lost_inside = 0
    counter%lost_outside = 0
    counter%tetr_pushings = 0
    counter%phi_0_mappings = 0

end subroutine set_local_counter_zero

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

subroutine initialize_poloidal_flux(poloidal_flux,num_particles)

    use boltzmann_types_mod, only: poloidal_flux_t
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: ntetr, tetra_grid

    integer, intent(in) :: num_particles
    integer :: i
    type(poloidal_flux_t), intent(out) :: poloidal_flux

    allocate(poloidal_flux%particle(num_particles))
    poloidal_flux%particle = 0
    poloidal_flux%max = 0
    poloidal_flux%min = tetra_physics(1)%Aphi1
    do i = 1, ntetr
        poloidal_flux%max = max(poloidal_flux%max,tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
                            & (verts([1,2,3],tetra_grid(i)%ind_knot(4))-verts([1,2,3],tetra_grid(i)%ind_knot(1)))))
        poloidal_flux%min = min(poloidal_flux%min,tetra_physics(i)%Aphi1)
    enddo

end subroutine initialize_poloidal_flux

subroutine calc_starting_conditions(b, v0, start_pos_pitch_mat, randcol)
    
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, ntetr
    use find_tetra_mod, only: find_tetra
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_physics_mod, only: coord_system
    use boltzmann_types_mod, only: boltzmann_input_t

    type(boltzmann_input_t), intent(in)                    :: b
    real(dp), intent(in)                                   :: v0
    real(dp), dimension(:,:), allocatable, intent(out)     :: start_pos_pitch_mat
    real(dp), dimension(:,:,:)                             :: randcol
    real(dp)                                               :: rand_scalar, vpar, vperp
    real(dp)                                               :: amin, cmin, cmax !amax is set globally
    real(dp), dimension(:), allocatable                    :: rand_vector
    real(dp), dimension(:,:), allocatable                  :: rand_matrix1, rand_matrix2
    integer                                                :: i
    real(dp), dimension(3)                                 :: x
    integer                                                :: ind_tetr_out,iface,seed_inp_unit
    real(dp)                                               :: amax, constant_part_of_weights
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

    allocate(start_pos_pitch_mat(5,b%num_particles))
    allocate(rand_vector(b%num_particles))
    allocate(rand_matrix1(3,b%num_particles))
    allocate(rand_matrix2(2,b%num_particles))

    start_pos_pitch_mat = 0

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

    constant_part_of_weights = b%density*(amax-amin)*(cmax-cmin)*2*pi

    !compute starting conditions
    if (b%boole_point_source) then
        if (grid_kind.eq.2) then
            start_pos_pitch_mat(1,:) = 160
            start_pos_pitch_mat(2,:) = 0
            start_pos_pitch_mat(3,:) = 70
        elseif (grid_kind.eq.4) then
            start_pos_pitch_mat(1,:) = 205
            start_pos_pitch_mat(2,:) = 0
            start_pos_pitch_mat(3,:) = 0
        endif
        if (coord_system.eq.2) print*, 'error: point source is only implemented for cylindrical coordinate system'
    else
        call RANDOM_NUMBER(rand_matrix1)
        start_pos_pitch_mat(ind_a,:) = amin + (amax - amin)*rand_matrix1(ind_a,:) !r in cylindrical, s in flux coordinates
        start_pos_pitch_mat(ind_b,:) = 2*pi*rand_matrix1(ind_b,:) !phi in cylindrical and flux coordinates
        start_pos_pitch_mat(ind_c,:) = cmin + (cmax - cmin)*rand_matrix1(ind_c,:) !z in cylindrical, theta in flux coordinates
        ! start_pos_pitch_mat(ind_a,:) = (/(214 + i*(216-214)/n_particles, i=1,b%num_particles)/)!r
        ! start_pos_pitch_mat(ind_b,:) = 0.0d0  !1d-1 !phi in cylindrical and flux coordinates
        ! start_pos_pitch_mat(ind_c,:) = 12d0 !z in cylindrical, theta in flux coordinates
    endif

    call RANDOM_NUMBER(rand_matrix2)
    !start_pos_pitch_mat(4,:) = 2*rand_matrix2(1,:)-1 !pitch parameter
    start_pos_pitch_mat(4,:) = 0.7d0 !pitch parameter

    if (b%boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start_pos_pitch_mat(5,:) = 5*b%energy_eV*rand_matrix2(2,:) !boltzmann energy distribution
        constant_part_of_weights = constant_part_of_weights*10/sqrt(pi)*b%energy_eV*ev2erg
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (b%boole_antithetic_variate) then
        start_pos_pitch_mat(:,1:b%num_particles:2) = start_pos_pitch_mat(:,2:b%num_particles:2)
        start_pos_pitch_mat(4,1:b%num_particles:2) = -start_pos_pitch_mat(4,2:b%num_particles:2)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    weights = constant_part_of_weights
    if (b%boole_refined_sqrt_g.eqv..false.) weights = constant_part_of_weights*start_pos_pitch_mat(ind_a,:)

    if (b%boole_precalc_collisions) then
        call RANDOM_NUMBER(randcol)
        !3.464102 = sqrt(12), this creates a random number with zero average and unit variance
        randcol(:,:,1:2:3) =  3.464102*(randcol(:,:,1:2:3)-.5)
    endif
    print*, 'calc_starting_conditions finished'
end subroutine calc_starting_conditions

subroutine prepare_collisions(c,v0,energy_eV)

    use boltzmann_types_mod, only: collisions_t
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_physics_mod, only: particle_mass,particle_charge
    use constants, only: ev2erg, echarge,amp
    use tetra_grid_settings_mod, only: grid_size
    use collis_ions, only: collis_init

    real(dp), intent(in) :: v0, energy_eV
    type(collisions_t) :: c
    integer :: n
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat
    integer :: aTe_unit, aTi_unit, ane_unit
    integer :: i, j
    real(dp) :: m0, z0

    n = 2 !number of background species
    allocate(c%dens_mat(n-1,ntetr))
    allocate(c%temp_mat(n,ntetr))
    allocate(c%vpar_mat(n,ntetr))
    allocate(c%efcolf_mat(n,ntetr))
    allocate(c%velrat_mat(n,ntetr))
    allocate(c%enrat_mat(n,ntetr))
    allocate(c%mass_num(n-1))
    allocate(c%charge_num(n-1))
    allocate(c%dens(n))
    allocate(c%temp(n))
    allocate(efcolf(n))
    allocate(velrat(n))
    allocate(enrat(n))
    c%mass_num = 0
    c%charge_num = 0
    c%mass_num(1) = 2
    !c%mass_num(2) = 3
    c%charge_num(1) = 1
    !c%charge_num(2) = 2
    c%vpar_mat = 0 !ask Sergei when this will be needed!!!
    m0 = particle_mass/amp
    z0 = particle_charge/echarge

    open(newunit = aTe_unit, file = 'background/Te_d.dat')
    read(aTe_unit,'(e16.9)') (c%temp_mat(2,i),i=1,ntetr/grid_size(2),3)
    close(aTe_unit)

    open(newunit = aTi_unit, file = 'background/Ti_d.dat')
    read(aTi_unit,'(e16.9)') (c%temp_mat(1,i),i=1,ntetr/grid_size(2),3)
    close(aTi_unit)

    open(newunit = ane_unit, file = 'background/ne_d.dat')
    read(ane_unit,'(e16.9)') (c%dens_mat(1,i),i=1,ntetr/grid_size(2),3)
    close(ane_unit)

    do i = 1,grid_size(2)-1 !copy data from first phi slice to all other phi slices
        c%temp_mat(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = c%temp_mat(:,1:ntetr/grid_size(2):3)
        c%dens_mat(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = c%dens_mat(:,1:ntetr/grid_size(2):3)
    enddo
    do i = 1,2 !copy data from first tetrahedron of each triangular prism to the two other ones
        c%temp_mat(:,1+i:ntetr:3) = c%temp_mat(:,1:ntetr:3)
        c%dens_mat(:,1+i:ntetr:3) = c%dens_mat(:,1:ntetr:3)
    enddo

    do i = 1, ntetr
        do j = 1,n
            !> if statement because electron density will be calculated in collis init
            if (j.lt.n) c%dens(j) = c%dens_mat(j,i)
            c%temp(j) = c%temp_mat(j,i)
        enddo
        call collis_init(m0,z0,c%mass_num,c%charge_num,c%dens,c%temp,energy_eV,v0,efcolf,velrat,enrat)
        c%efcolf_mat(:,i) = efcolf
        c%velrat_mat(:,i) = velrat
        c%enrat_mat(:,i) = enrat
    enddo
end subroutine prepare_collisions

end module boltzmann_mod