!
!===============================================================================
! alpha_lifetime_gorilla_mod.f90
!
! Compatibility client for the reusable marker-transport library
! (SRC/transport).  This applet traces alpha-particle guiding-centre orbits and
! records their confinement time.  The orbit integration itself is delegated
! to GORILLA's characteristic-step stepper through the library's generic
! advance_ensemble API; the module only owns the applet input parsing, the
! start-state sampling and the output formatting.
!
! The i_integrator_type == 1 (GORILLA pusher) path is migrated to the library.
! The type 0 (no orbit) and type 2 (canonical) paths are preserved as-is.
!===============================================================================
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
            integer :: inp_unit
!
            open(newunit=inp_unit, file='alpha_lifetime.inp', status='unknown')
            read(inp_unit,nml=alpha_lifetimenml)
            close(inp_unit)

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
            integer :: j,n,seed_unit
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
                    open(newunit=seed_unit,file='seed.inp')
                    write(seed_unit,*) n
                    write(seed_unit,*) seed
                    close(seed_unit)
                case(2) !Load seed
                    open(newunit = seed_unit, file='seed.inp', status='old',action = 'read')
                    read(seed_unit,*) n
                    allocate(seed(n))
                    read(seed_unit,*) seed
                    close(seed_unit)
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
!       Library stepper wrapper: adapts GORILLA's orbit_timestep_gorilla
!       characteristic-step contract to the library step_signature_t.
!
        subroutine alpha_gorilla_step(x,vpar,vperp,dt,boole_initialized,cell,iface,t_remain,boole_lost)
!
            use orbit_timestep_gorilla_mod, only: orbit_timestep_gorilla
            use, intrinsic :: iso_fortran_env, only: dp => real64
!
            implicit none
!
            real(dp), dimension(3), intent(inout) :: x
            real(dp), intent(inout) :: vpar,vperp
            real(dp), intent(in)    :: dt
            logical,  intent(inout) :: boole_initialized
            integer,  intent(inout) :: cell,iface
            real(dp), intent(out)   :: t_remain
            logical,  intent(out)   :: boole_lost
!
            t_remain = dt
            call orbit_timestep_gorilla(x,vpar,vperp,dt,boole_initialized,cell,iface,t_remain)
            boole_lost = (cell.eq.-1)
!
        end subroutine alpha_gorilla_step
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine calc_alpha_lifetime_gorilla()
!
            use orbit_timestep_gorilla_mod, only: initialize_gorilla
            use constants, only: ev2erg, pi
            use tetra_physics_mod, only: particle_mass,cm_over_e,mag_axis_R0
            use fluxtv_mod, only: load_flux_tube_volume,pos_fluxtv_mat
            use parmot_mod, only: rmu,ro0
            use velo_mod, only: isw_field_type
            use supporting_functions_mod, only: theta_sym_flux2theta_vmec
            use tetra_grid_settings_mod, only: n_field_periods
            use tetra_grid_mod, only: ntetr
            use binsrc_mod, only: binsrc
            use sub_alpha_lifetime_can_mod, only: orbit_timestep_can
            use marker_transport_mod, only: ensemble_t,advance_ensemble,noop_process, &
                & transport_config_t,termination_event_t,TERM_LOST
            use conservative_tallies_mod, only: tallies_t
            use, intrinsic :: iso_fortran_env, only: int64, dp => real64
!
            implicit none
!
            double precision :: vmod,pitchpar,vpar,vperp,t_remain,t_confined,tau_out_can
            integer :: kpart,i,n,ind_tetr,iface,n_lost_particles,ierr
            integer :: n_start, n_end, i_part, alpha_unit
            double precision, dimension(3) :: x_rand_beg,x
            double precision, dimension(:), allocatable :: xi
            double precision, dimension(:), allocatable :: t_conf
            logical :: boole_initialized,boole_particle_lost
            double precision :: dtau, dphi,dtaumin
            double precision, dimension(5) :: z
!
            !Transport-library state
            type(ensemble_t)              :: ens
            type(tallies_t)               :: tallies
            type(transport_config_t)      :: cfg
            type(termination_event_t), allocatable :: evts(:)
            real(dp), dimension(:), allocatable :: cell_volume
            real(dp) :: init_weight, init_momentum, init_energy
!
            !Load input for alpha lifetime computation
            call load_alpha_lifetime_inp()
!
            allocate(xi(n_particles))
            allocate(t_conf(n_particles))
            t_conf = 0.d0
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
            open(newunit=alpha_unit,file=filename_alpha_lifetime)
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
            n_lost_particles = 0
!
            select case(i_integrator_type)
!
!===============================================================================
! Migrated path (i_integrator_type == 1): GORILLA pusher through the library.
!===============================================================================
                case(1)
!
                    !--- Build marker ensemble (start states sampled as before) ---
                    call ens%init(n_particles, 0_int64)
!
                    do n = n_start,n_end
                        !Find random start indices that are distributed proportionally to the flux tube volume
                        call binsrc(pos_fluxtv_mat(:,4),1,size(pos_fluxtv_mat(:,4)),xi(n),i)
!
                        x_rand_beg = pos_fluxtv_mat(i,1:3)
!
                        ens%markers(n)%x = x_rand_beg
!
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
                        vpar = pitchpar * vmod
                        vperp = sqrt(vmod**2-vpar**2)
!
                        ens%markers(n)%vpar  = vpar
                        ens%markers(n)%vperp = vperp
                        ens%markers(n)%pitch = pitchpar
                        ens%markers(n)%weight = 1.d0
                        ens%markers(n)%cell   = -1
                        ens%markers(n)%iface  = -1
                        ens%markers(n)%boole_initialized = .false.
                        ens%markers(n)%active = .true.
                    enddo
!
                    !--- Diagnostic tallies (cell/surface + conservation ledger) ---
                    call tallies%cells%init(ntetr)
                    call tallies%surfaces%init(ntetr)
                    allocate(cell_volume(ntetr))
                    cell_volume = 1.d0
                    call tallies%cells%norm%init(1.d0, time_step, n_particles, cell_volume)
                    call tallies%surfaces%norm%init(1.d0, time_step, n_particles)
                    call ens%totals(init_weight, init_momentum, init_energy)
                    call tallies%ledger%begin_from_totals(init_weight, init_momentum, init_energy)
!
!                   --- Configure and run the library transport ---
                    allocate(evts(n_particles))
                    cfg%n_steps = 1
                    cfg%dt      = time_step
                    cfg%boole_record_cell_tallies      = .true.
                    cfg%boole_record_surface_crossings = .true.
!
                    call advance_ensemble(ens, alpha_gorilla_step, noop_process, &
                                          tallies, cfg, evts)
!
                    !--- Collect confinement times and lost count ---
                    do n = 1,n_particles
                        t_conf(n) = ens%markers(n)%time
                        if(evts(n)%occurred .and. evts(n)%event_type.eq.TERM_LOST) then
                            n_lost_particles = n_lost_particles + 1
                        endif
                    enddo
!
                    !--- Write confinement times (deterministic marker order) ---
                    do n = 1,n_particles
                        write(alpha_unit,*) t_conf(n)
                    enddo
!
                case(0)
                    !No orbit computation: all particles confined for the full step
                    do n = 1,n_particles
                        t_conf(n) = time_step
                        write(alpha_unit,*) t_conf(n)
                    enddo
!
                case(2)
!
                    !Legacy canonical (direct) integrator path, preserved as-is.
                    do n = n_start,n_end
!
                        !Find random start indices that are distributed proportionally to the flux tube volume
                        call binsrc(pos_fluxtv_mat(:,4),1,size(pos_fluxtv_mat(:,4)),xi(n),i)
!
                        x_rand_beg = pos_fluxtv_mat(i,1:3)
!
                        z(1) = x_rand_beg(1)
                        z(3) = x_rand_beg(3)
                        z(2) = theta_sym_flux2theta_vmec(z(1),x_rand_beg(2),z(3))  !Transform theta_symflux to theta_vmec
                        z(4) = 1.d0
!
                        !Find random pitch parameter
                        if(boole_random_precalc) then
                            pitchpar = rd_start_pitchpar(n)
                        else
                            !$omp critical
                                call random_number(pitchpar)
                            !$omp end critical
                        endif
                        pitchpar = 2.d0*pitchpar - 1.d0
                        z(5) = pitchpar
!
                        ierr = 0
                        call orbit_timestep_can(z,dtau,dtaumin,ierr,tau_out_can)
!
                        if(ierr.eq.1) then
                            n_lost_particles = n_lost_particles + 1
                        endif
!
                        t_confined = -1.d0*tau_out_can/vmod
!
                    enddo
!
            end select
!
            close(alpha_unit)
print *, 'Number of lost particles',n_lost_particles
!
            !Deallocate random numbers
            if(boole_random_precalc) then
                call dealloc_rand_numbers_alpha_lifetime()
            endif
!
        end subroutine calc_alpha_lifetime_gorilla
!
    end module alpha_lifetime_gorilla_mod
