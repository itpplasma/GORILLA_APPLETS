!
    module alpha_lifetime_gorilla_mod
!
    implicit none
!
        private
!
        double precision, dimension(:), allocatable     :: rd_start_position
        double precision, dimension(:), allocatable     :: rd_start_pitchpar
!
        double precision :: time_step,energy_eV
        integer :: n_particles, i_integrator_type, seed_option
        logical :: boole_random_precalc
        character(1024) :: filename_alpha_lifetime
!
        !Namelist for Alpha Lifetime input
        NAMELIST /alpha_lifetimenml/ time_step,energy_eV,n_particles, i_integrator_type, &
            & boole_random_precalc,filename_alpha_lifetime,seed_option
!
        public :: calc_alpha_lifetime_gorilla
!
    contains
!
        subroutine load_alpha_lifetime_inp()
!
            open(unit=71, file='alpha_lifetime.inp', status='unknown')
            read(71,nml=alpha_lifetimenml)
            close(71)

            print *,'GORILLA: Loaded input data from alpha_lifetime.inp'
!
        end subroutine load_alpha_lifetime_inp
!
        subroutine calc_rand_numbers_alpha_lifetime(n_particles)
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

        end subroutine calc_rand_numbers_alpha_lifetime
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine dealloc_rand_numbers_alpha_lifetime()
!
            deallocate(rd_start_position)
            deallocate(rd_start_pitchpar)
!
        end subroutine dealloc_rand_numbers_alpha_lifetime
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine calc_alpha_lifetime_gorilla()
!
            use orbit_timestep_gorilla_mod, only: initialize_gorilla,orbit_timestep_gorilla
            use constants, only: ev2erg, pi
            use tetra_physics_mod, only: particle_mass,cm_over_e,mag_axis_R0
            use fluxtv_mod, only: load_flux_tube_volume,pos_fluxtv_mat
            !use omp_lib, only: omp_get_thread_num
            use parmot_mod, only: rmu,ro0
            use velo_mod, only: isw_field_type
            use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
            use tetra_grid_settings_mod, only: n_field_periods
!
            implicit none
!
            double precision :: vmod,pitchpar,vpar,vperp,t_remain,t_confined,tau_out_can
            integer :: kpart,i,n,ind_tetr,iface,n_lost_particles,ierr
            integer :: n_start, n_end, i_part
            double precision, dimension(3) :: x_rand_beg,x
            double precision, dimension(:), allocatable :: xi
            logical :: boole_initialized,boole_particle_lost
            double precision :: dtau, dphi,dtaumin
            double precision, dimension(5) :: z
!
            !Load input for alpha lifetime computation
            call load_alpha_lifetime_inp()
!
            allocate(xi(n_particles))
!
            n_start = 1
            n_end = n_particles
!
            !Initialize GORILLA
            call initialize_gorilla()
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
            !Open file for writing alpha_life_time_gorilla
            open(99,file=filename_alpha_lifetime)
!

            !Precompute random numbers
            if(boole_random_precalc) then
                call calc_rand_numbers_alpha_lifetime(n_particles)
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
            !$OMP& dtau,dtaumin,rd_start_position,rd_start_pitchpar,boole_random_precalc,n_start,n_end) &
            !$OMP& PRIVATE(n,boole_particle_lost,i,x_rand_beg,x,pitchpar,vpar,vperp,boole_initialized, &
            !$OMP& ind_tetr,iface,t_remain,t_confined,z,ierr,tau_out_can)
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
!
                !Orbit integration
                select case(i_integrator_type)
!
                    case(0)
                        !No orbit computation
!
                    case(1)
                        boole_initialized = .false.
                        call orbit_timestep_gorilla(x,vpar,vperp,time_step,boole_initialized,ind_tetr,iface,t_remain)
!
                        !Confinement time of alpha particle
                        t_confined = time_step - t_remain
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
                            !write(99,*) n, boole_particle_lost , x_rand_beg ,pitchpar,x(1),t_confined
                            write(99,*) t_confined
                        !$omp end critical
!
                    case(2)
                        ierr = 0
                        call orbit_timestep_can(z,dtau,dtaumin,ierr,tau_out_can)
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
                        !Write results in file
!                        !$omp critical
!                            write(99,*) n, boole_particle_lost , x_rand_beg ,pitchpar,z(1),t_confined
!                        !$omp end critical
!
                end select
!
            enddo !n
            !$OMP END DO
            !$OMP END PARALLEL
!
!            close(99)
print *, 'Number of lost particles',n_lost_particles
!open(99,file='confined_fraction.dat')
!write(99,*) 1.d0-dble(n_lost_particles)/dble(n_particles)
!
            !Deallocate random numbers
            if(boole_random_precalc) then
                call dealloc_rand_numbers_alpha_lifetime()
            endif
!
        end subroutine calc_alpha_lifetime_gorilla
!
    end module alpha_lifetime_gorilla_mod
