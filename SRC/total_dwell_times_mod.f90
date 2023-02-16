module total_dwell_times_mod
    implicit none

    private
    !
            double precision, dimension(:), allocatable     :: rd_start_position
            double precision, dimension(:), allocatable     :: rd_start_pitchpar
            

    double precision :: time_step,energy_eV
    integer :: n_particles, i_integrator_type, seed_option
    logical :: boole_random_precalc
    character(1024) :: filename_total_dwell_times, filename_starting_conditions
    !
    !Namelist for total dwell times input
    NAMELIST /total_dwell_timesnml/ time_step,energy_eV,n_particles, i_integrator_type, &
    & boole_random_precalc,filename_total_dwell_times,seed_option,filename_starting_conditions

    public :: calc_total_dwell_times
    
contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine load_total_dwell_times_inp()
    !
    open(unit=71, file='total_dwell_times.inp', status='unknown')
    read(71,nml=total_dwell_timesnml)
    close(71)
    
    print *,'GORILLA: Loaded input data from total_dwell_times.inp'
end subroutine load_total_dwell_times_inp

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine calc_rand_numbers_total_dwell_times(n_particles)
    !
    implicit none
!
    integer, intent(in)     :: n_particles
    integer,dimension(:), allocatable :: seed
    double precision, dimension(:), allocatable :: rd_seed
    integer :: j,n
!
    allocate(rd_start_position(n_particles))
    allocate(rd_start_pitchpar(n_particles))
!
!           !seed_option (Input file: gorilla_applets.inp)
    ! 1 ... produce seed, 2 ... load seed
    select case(seed_option)
        case(1) !Produce seed
            !Allocate seed
            n = 0
            call random_seed(size=n)
            allocate(seed(n))
!
            allocate(rd_seed(n))
            !$omp critical
                call random_number(rd_seed)
            !$omp end critical
            seed = int(rd_seed*10.d0)
            deallocate(rd_seed)
!
            open(85,file='seed.inp')
            write(85,*) n
            write(85,*) seed
            close(85)
        case(2) !Load seed
            open(unit = 85, file='seed.inp', status='old',action = 'read')
            read(85,*) n
            allocate(seed(n))
            read(85,*) seed
            close(85)
    end select
!
    CALL RANDOM_SEED (PUT=seed)
!
    deallocate(seed)
!
    !$omp critical
        call random_number(rd_start_position)
        call random_number(rd_start_pitchpar)
    !$omp end critical
    
end subroutine calc_rand_numbers_total_dwell_times

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine dealloc_rand_numbers_total_dwell_times()
    !
    deallocate(rd_start_position)
    deallocate(rd_start_pitchpar)
    !
end subroutine dealloc_rand_numbers_total_dwell_times

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine calc_total_dwell_times

    use orbit_timestep_gorilla_mod, only: initialize_gorilla !,orbit_timestep_gorilla
    use constants, only: ev2erg, pi
    use tetra_physics_mod, only: particle_mass,particle_charge,cm_over_e,mag_axis_R0
    use fluxtv_mod, only: load_flux_tube_volume,pos_fluxtv_mat
    !use omp_lib, only: omp_get_thread_num
    use parmot_mod, only: rmu,ro0
    use velo_mod, only: isw_field_type
    use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
    use tetra_grid_settings_mod, only: n_field_periods
    use tetra_grid_mod, only: ntetr

    implicit none

    double precision, dimension(:), allocatable :: total_dwell_times, single_particle_dwell_times, densities
    double precision, dimension(:,:), allocatable :: total_currents, single_particle_currents, fluid_velocities, start_pos_pitch_mat
    double precision :: vmod,pitchpar,vpar,vperp,t_remain,t_confined,tau_out_can
    integer :: kpart,i,n,l,m,ind_tetr,iface,n_lost_particles,ierr
    integer :: n_start, n_end, i_part
    double precision, dimension(3) :: x_rand_beg,x
    double precision, dimension(:), allocatable :: xi
    logical :: boole_initialized,boole_particle_lost
    double precision :: dtau, dphi,dtaumin
    double precision, dimension(5) :: z

    !Load input for total dwell times computation
    call load_total_dwell_times_inp()
    !
    allocate(xi(n_particles))
!
    n_start = 1
    n_end = n_particles
!
    !Initialize GORILLA
    call initialize_gorilla()

    allocate(total_dwell_times(1:ntetr))
    allocate(single_particle_dwell_times(1:ntetr))
    allocate(densities(1:ntetr))
    allocate(total_currents(1:3,1:ntetr))
    allocate(single_particle_currents(1:3,1:ntetr))
    allocate(fluid_velocities(1:3,1:ntetr))
    total_dwell_times = 0
    total_currents = 0
!
    !Load fluxtube volume for a starting position (File to be chosen in gorilla_applets.inp)
    call load_flux_tube_volume()
!
    !Compute velocity module from kinetic energy dependent on particle species
    vmod=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
    !
    !------------------------------------------------------------------------------------------------------------!
    !------------------------------------ Initialization of direct integrator -----------------------------------!
    !
                if(i_integrator_type.eq.2) then
    !
                    !inverse relativistic temperature
                    rmu=1d8
    !
                    !normalized larmor radius
                    ro0 = vmod*cm_over_e
    !
                    isw_field_type=1
    
                    !normalized slowing down time:
                    dtau = -1.d0*time_step*vmod
    !
                    !field line integration step step over phi (to check chamber wall crossing)
                    dphi=2.d0*pi/dble(100)
    !
                    !orbit integration time step (to check chamber wall crossing)
                    dtaumin=dphi*mag_axis_R0
    !
                endif
    !
    !------------------------------------------------------------------------------------------------------------!
    !
                !Open file for writing total dwell times
                !open(99,file=filename_total_dwell_times)
    !
    
                !Precompute random numbers
                if(boole_random_precalc) then
                    call calc_rand_numbers_total_dwell_times(n_particles)
                endif
    !
                !Create random numbers for sampling starting points on flux surface
                !Flux tube volume for certain flux-sufrace needs to be computed in previous compilation step
                if(boole_random_precalc) then
                    xi = rd_start_position
                else
                    !$omp critical
                        call random_number(xi)
                    !$omp end critical
                endif
    !
                kpart = 0
                n_lost_particles = 0
    !
                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP& SHARED(n_particles,pos_fluxtv_mat,xi,n_lost_particles,kpart,vmod,time_step,i_integrator_type, &
                !$OMP& dtau,dtaumin,rd_start_position,rd_start_pitchpar,boole_random_precalc,n_start,n_end, total_dwell_times, &
                !$OMP& total_currents, start_pos_pitch_mat) &
                !$OMP& PRIVATE(n,boole_particle_lost,i,x_rand_beg,x,pitchpar,vpar,vperp,boole_initialized, &
                !$OMP& ind_tetr,iface,t_remain,t_confined,z,ierr,tau_out_can, single_particle_dwell_times, &
                !$OMP& single_particle_currents)
                !$OMP DO
    !
                !Loop over particles
                do n = n_start,n_end !1,n_particles
    !
                    !Counter for particles
                    !$omp critical
                    kpart = kpart+1
                    boole_particle_lost = .false.
    print *, kpart, ' / ', n_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
                    !$omp end critical
    
                    !Find random start indices that are distributed proportionally to the flux tube volume
                    call binsrc(pos_fluxtv_mat(:,4),1,shape(pos_fluxtv_mat(:,4)),xi(n),i)
    !
                    x_rand_beg = pos_fluxtv_mat(i,1:3)
    !
                    select case(i_integrator_type)
                        case(1,0)
                            x = x_rand_beg
                        case(2)
                            z(1) = x_rand_beg(1)
                            z(3) = x_rand_beg(3)
                            z(2) = theta_sym_flux2theta_vmec(z(1),x_rand_beg(2),z(3))  !Transform theta_symflux to theta_vmec
                            z(4) = 1.d0
                    end select
    
                    !Find random pitch parameter
                    if(boole_random_precalc) then
                        pitchpar = rd_start_pitchpar(n)
                    else
                        !$omp critical
                            call random_number(pitchpar)
                        !$omp end critical
                    endif
                    pitchpar = 2.d0*pitchpar - 1.d0
    !
                    !Particle velocities in accordance with integrator
                    select case(i_integrator_type)
                        case(1,0)
                            vpar = pitchpar * vmod
                            vperp = sqrt(vmod**2-vpar**2)
                        case(2)
                            z(5) = pitchpar
                    end select
                    
                    !Set single_particle_dwell_times and single_particle_currents back to zero for every loop iteration
                    single_particle_dwell_times = 0
                    single_particle_currents = 0
    !
                    !Orbit integration
                    select case(i_integrator_type)
    !
                        case(0)
                            !No orbit computation
    !
                        case(1)
                            boole_initialized = .false.
                            call orbit_timestep_gorilla_dwell_times(x,vpar,vperp,time_step,boole_initialized, &
                            & ind_tetr,iface,t_remain, single_particle_dwell_times, single_particle_currents)
    !
                            !Confinement time of particles
                            t_confined = time_step - t_remain

                            !dwell times of particles

                            !$omp critical
                            total_dwell_times = total_dwell_times + single_particle_dwell_times
                            total_currents = total_currents + single_particle_currents
                            !$omp end critical
    !
                            !Lost particle handling
                            if(ind_tetr.eq.-1) then
                                !$omp critical
                                    n_lost_particles = n_lost_particles + 1
                                    boole_particle_lost = .true.
                                !$omp end critical
                            endif
    !
                            !Write results in file
                            !$omp critical
                                write(99,*) n, boole_particle_lost , x_rand_beg ,pitchpar,x(1),t_confined
                            !$omp end critical
    !
                        case(2)
                            ierr = 0
                            call orbit_timestep_can_dwell_times(z,dtau,dtaumin,ierr,tau_out_can)
    !
                            if(ierr.eq.1) then
                                !$omp critical
                                    n_lost_particles = n_lost_particles + 1
                                    boole_particle_lost = .true.
                                !$omp end critical
                            endif
    !
                            t_confined = -1.d0*tau_out_can/vmod
    !
    !                        !Write results in file
                            !$omp critical
                                write(99,*) n, boole_particle_lost , x_rand_beg ,pitchpar,z(1),t_confined
                            !$omp end critical
    !
                    end select
    !
                enddo !n
                !$OMP END DO
                !$OMP END PARALLEL
    !
    !            close(99)
                
                densities = total_dwell_times / time_step !compute densities
                do l = 1,ntetr !compute fluid velocity
                        if (densities(l) .EQ. 0) then
                            fluid_velocities(:,l) = 0
                        else
                            do m = 1,3
                                fluid_velocities(m,l) = total_currents(m,l)/(densities(l)*particle_charge)
                            enddo
                        endif
                enddo


    print *, 'Number of lost particles',n_lost_particles
    open(99,file='confined_fraction.dat')
    write(99,*) 1.d0-dble(n_lost_particles)/dble(n_particles)
    !
                !Deallocate random numbers
                if(boole_random_precalc) then
                    call dealloc_rand_numbers_total_dwell_times()
                endif

    
    open(50, file = filename_total_dwell_times)
    write(50,'(ES15.3E4)') total_dwell_times
    close(50)
    open(51, file = 'densities.dat')
    write(51,'(ES15.3E4)') densities
    close(51)
    open(52, file = 'total_currents.dat')
    write(52,'(3ES15.3E4)') total_currents
    close(52)
    open(53, file = 'fluid_velocities.dat')
    write(53,'(3ES15.3E4)') fluid_velocities
    close(53)

                deallocate(total_dwell_times, single_particle_dwell_times, total_currents, &
                & single_particle_currents, fluid_velocities, densities)

PRINT*, 'particle mass = ', particle_mass
PRINT*, 'large radius = ', mag_axis_R0
PRINT*, 'parallel velocity = ', vpar
PRINT*, 'absolute value of velocity = ', vmod
PRINT*, 'perpendicular velocity = ', vperp
PRINT*, 'pitch par =', pitchpar
PRINT*, 'particle charge = ', particle_charge
    !
end subroutine calc_total_dwell_times

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    

subroutine orbit_timestep_gorilla_dwell_times(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, &
    & t_remain_out,single_particle_dwell_times, single_particle_currents)
    !
                use pusher_tetra_rk_mod, only: find_tetra,pusher_tetra_rk,initialize_const_motion_rk
                use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
                use tetra_physics_poly_precomp_mod , only: make_precomp_poly_perpinv, initialize_boole_precomp_poly_perpinv, &
                    & alloc_precomp_poly_perpinv
                use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass
                !use tetra_grid_mod, only: ntetr
                use gorilla_settings_mod, only: ipusher, poly_order
                use orbit_timestep_gorilla_mod, only: check_coordinate_domain
                use supporting_functions_mod, only: bmod_func, vperp_func
    !
                implicit none
    !
                double precision, dimension(3), intent(inout)   :: x
                double precision, intent(inout)                 :: vpar,vperp
                double precision, intent(in)                    :: t_step
                logical, intent(inout)                          :: boole_initialized
                integer, intent(inout)                          :: ind_tetr,iface
                double precision, intent(out), optional         :: t_remain_out
                double precision, dimension(3)                  :: z_save, x_save
                double precision                                :: vperp2,t_remain,t_pass,vpar_save
                logical                                         :: boole_t_finished
                integer                                         :: ind_tetr_save,iper,k
                double precision                                :: perpinv,perpinv2
                double precision, dimension(:), intent(inout)   :: single_particle_dwell_times
                double precision, dimension(:,:), intent(inout) :: single_particle_currents
    !
                !If orbit_timestep is called for the first time without grid position
                if(.not.boole_initialized) then
    !
                    !Check coordinate domain (optionally perform modulo operation)
                    call check_coordinate_domain(x)
    !
                    !Find tetrahedron index and face index for position x
                    call find_tetra(x,vpar,vperp,ind_tetr,iface)
    !               
                    !If particle doesn't lie inside any tetrahedron
                    if(ind_tetr.eq.-1) then
                        return
                    endif
    !
                    boole_initialized = .true.
                endif
    !           
                !Exit the subroutine after initialization, if time step equals zero
                if(t_step.eq.0.d0) return
    !
                !Squared perpendicular velocity
                vperp2 = vperp**2
    !
                !Compute relative particle position
                z_save = x-tetra_physics(ind_tetr)%x1
    !
                !Compute perpendicular invariant of particle
                perpinv=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)
                perpinv2 = perpinv**2
    !               
                !Initialize constants of motion in particle-private module
                select case(ipusher)
                    case(1)
                        call initialize_const_motion_rk(perpinv,perpinv2)
                    case(2)
                        call initialize_const_motion_poly(perpinv,perpinv2)
                end select        
    !
    !             !NOT FULLY IMPLEMENTED YET: Precompute quatities dependent on perpinv
    !             call alloc_precomp_poly_perpinv(1,ntetr)
    !             call initialize_boole_precomp_poly_perpinv()
    !             call make_precomp_poly_perpinv(perpinv,perpinv2)
    !
                !Integrate particle orbit for given time step
                t_remain = t_step
    !
                !Logical for handling time integration
                boole_t_finished = .false.
    !
                !Loop for tetrahedron pushings until t_step is reached
                do
    !
                    !Domain Boundary
                    if(ind_tetr.eq.-1) then
                        print *, 'WARNING: Particle lost.'
                        if( present(t_remain_out)) then
                            t_remain_out = t_remain
                        endif
                        exit
                    endif
    !                
                    !Save the tetrahedron index for computation of vperp in the last step
                    ind_tetr_save = ind_tetr
    !
                    !Save vpar for the computation of the parallel adiabatic invariant
                    vpar_save = vpar
    !  
                    !t_remain (in) ... remaining time until t_step is finished
                    !t_pass (out) ... time to pass the tetrahdron

                    !Save x for the computation of the current
                    x_save = x
    !
                    !Calculate trajectory
                    select case(ipusher)
                        case(1)
                            call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper)
                        case(2)
                            call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_remain,&
                                                                & t_pass,boole_t_finished,iper)
                    end select
    !
                    t_remain = t_remain - t_pass

                    single_particle_dwell_times(ind_tetr_save) = single_particle_dwell_times(ind_tetr_save) + t_pass
                    single_particle_currents(:,ind_tetr_save) = particle_charge*(x-x_save)/t_step ! /t_pass vor the speed but *t_pass to calculate the time averge together with /t_step
    !
                    !Orbit stops within cell, because "flight"-time t_step has finished
                    if(boole_t_finished) then
                        if( present(t_remain_out)) then
                            t_remain_out = t_remain
                        endif
                        exit
                    endif
    !
                enddo !Loop for tetrahedron pushings
!
                !Compute vperp from position
                vperp = vperp_func(z_save,perpinv,ind_tetr_save)
    !            
    !             !NOT FULLY IMPLEMENTED YET: Deallocate precomputed quantities dependent on perpinv
    !             call alloc_precomp_poly_perpinv(2,ntetr)
    !         
end subroutine orbit_timestep_gorilla_dwell_times


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine orbit_timestep_can_dwell_times(z,dtau,dtaumin,ierr,tau_out)
    !
          use odeint_mod, only: odeint_allroutines
          use runge_kutta_mod, only: runge_kutta_allroutines
          use gorilla_settings_mod, only: rel_err_ode45
    !
          implicit none
    !
          integer, parameter          :: ndim=5, nstepmax=10000000
          double precision            :: relerr
    !
          integer :: ierr,j
          double precision :: dtau,dtaumin,phi,tau1,tau2,tau_out
    !
          double precision, dimension(2)    :: y
          double precision, dimension(ndim) :: z
    !
          external velo_can
    !
          !Use relative error definition from input file
          relerr = rel_err_ode45
    !
          if(dtaumin*nstepmax.le.abs(dtau)) then
            ierr=2
            print *,'orbit_timestep: number of steps exceeds nstepmax'
            return
          endif
    !
          ierr=0
          y(1)=z(1)
          y(2)=z(2)
          phi=z(3)
    !
          call chamb_can(y,phi,ierr)
    !
          tau_out = 0.d0
    !
          if(ierr.ne.0) return
          tau1=0.d0
          tau2=sign(dtaumin,dtau)      
    !
          do while(abs(tau2).lt.abs(dtau))
    !
    !print *, '(1) before RK4', z(1)
            call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
    
    !        call runge_kutta_allroutines(4,z,ndim,tau1,tau2,velo_can)
    !print *, '(1) after RK4', z(1)
    !
            y(1)=z(1)
            y(2)=z(2)
            phi=z(3)
    !
            call chamb_can(y,phi,ierr)
    !
            tau_out = tau2
    !
            if(ierr.ne.0) return
            tau1=tau2
            tau2=tau2+sign(dtaumin,dtau)
          enddo
    !
          tau2=dtau
    !
    !print *, '(2) before RK4', z
          call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
    
    !      call runge_kutta_allroutines(4,z,ndim,tau1,tau2,velo_can)
    !print *, '(2) after RK4', z
    !
          y(1)=z(1)
          y(2)=z(2)
          phi=z(3)
    !
          call chamb_can(y,phi,ierr)
    !
          tau_out = tau2
    !
          if(ierr.ne.0) return
    !
          end subroutine orbit_timestep_can_dwell_times

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module total_dwell_times_mod


