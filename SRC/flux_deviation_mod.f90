!
    module flux_deviation_mod
!
        implicit none
!
        private
!
        double precision, dimension(:), allocatable     :: rd_start_position
        double precision, dimension(:), allocatable     :: rd_start_pitchpar
        double precision, dimension(:,:), allocatable   :: rd_xi_collision
!
        public :: calc_flux_deviation,accumulate_valid_moments,lorentz_pitch_step
!
    contains
!
        pure function lorentz_pitch_step(pitch,eps_collision,xi_sign) result(pitch_new)
!
            implicit none
!
            double precision, intent(in) :: pitch,eps_collision,xi_sign
            double precision :: pitch_new
!
            pitch_new = pitch*(1.d0-eps_collision) &
                        & + xi_sign*sqrt(max(0.d0,(1.d0-pitch**2)*eps_collision))
!
        end function lorentz_pitch_step
!
        subroutine accumulate_valid_moments(delta_s_squared,valid_steps,sum_delta_s_squared, &
                                            & sum_delta_s_fourth,n_contributors)
!
            implicit none
!
            double precision, dimension(:), intent(in) :: delta_s_squared
            integer(kind=8), intent(in) :: valid_steps
            double precision, dimension(:), intent(inout) :: sum_delta_s_squared,sum_delta_s_fourth
            integer, dimension(:), intent(inout) :: n_contributors
!
            if(valid_steps <= 0) return
            sum_delta_s_squared(1:valid_steps) = sum_delta_s_squared(1:valid_steps) &
                                                  & + delta_s_squared(1:valid_steps)
            sum_delta_s_fourth(1:valid_steps) = sum_delta_s_fourth(1:valid_steps) &
                                                 & + delta_s_squared(1:valid_steps)**2
            n_contributors(1:valid_steps) = n_contributors(1:valid_steps)+1
!
        end subroutine accumulate_valid_moments
!
        subroutine calc_rand_numbers(n_particles,n_time_steps)
!
            use mono_energetic_transp_coef_settings_mod, only: seed_option,random_seed_filename
!
            implicit none
!
            integer, intent(in)     :: n_particles
            integer(kind=8),intent(in) :: n_time_steps
            integer,dimension(:), allocatable :: seed
            double precision, dimension(:), allocatable :: rd_seed
            integer :: j,n,seed_unit
!
            allocate(rd_start_position(n_particles))
            allocate(rd_start_pitchpar(n_particles))
            allocate(rd_xi_collision(n_particles,n_time_steps))
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
                    open(newunit=seed_unit,file=trim(random_seed_filename),status='replace')
                    write(seed_unit,*) n
                    write(seed_unit,*) seed
                    close(seed_unit)
                case(2) !Load seed
                    open(newunit=seed_unit,file=trim(random_seed_filename),status='old',action='read')
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
                call random_number(rd_xi_collision)
            !$omp end critical

        end subroutine calc_rand_numbers
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine dealloc_rand_numbers()

            deallocate(rd_start_position)
            deallocate(rd_start_pitchpar)
            deallocate(rd_xi_collision)

        end subroutine
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine calc_flux_deviation(n_particles,n_time_steps,t_step,vmod,boole_collisions,boole_random_precalc,coll_freq, &
                                    & file_id_transp_diff_coef,file_id_psi2,file_id_std_psi2,i_integrator_type,idiffcoef_output, &
                                    & nu_star,file_id_loss_events,file_id_loss_summary,file_id_particle_histories)
!
            use omp_lib, only: omp_get_thread_num
            use fluxtv_mod, only: pos_fluxtv_mat
            use orbit_timestep_gorilla_mod, only: orbit_timestep_gorilla
            use tetra_physics_mod, only: coord_system,psi_rphiz,vector_potential_sthetaphi_vmec,cm_over_e,mag_axis_R0, &
                                       & vector_potential_sthetaphi
            use tetra_grid_settings_mod, only: grid_kind,grid_size
            use magdata_in_symfluxcoor_mod, only: psitor_max
            use constants, only: pi
            use gorilla_diag_mod, only: transp_coef_trigger
            use parmot_mod, only: rmu,ro0
            use velo_mod, only: isw_field_type
            use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
            use mono_energetic_transp_coef_settings_mod, only: lambda_start,boole_psi_mat,transport_backend
            use gorilla_settings_mod, only: boole_helical_pert, helical_pert_eps_Aphi, helical_pert_m_fourier, &
                                          & helical_pert_n_fourier
            use binsrc_mod, only: binsrc
            use llsq_mod, only: llsq
            use sub_alpha_lifetime_can_mod, only: orbit_timestep_can
            use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
!
            implicit none
!
            integer,intent(in) :: n_particles,file_id_psi2,file_id_std_psi2,i_integrator_type,file_id_transp_diff_coef
            integer,intent(in) :: idiffcoef_output
            integer,intent(in),optional :: file_id_loss_events,file_id_loss_summary,file_id_particle_histories
            integer(kind=8),intent(in) :: n_time_steps
            double precision, intent(in) :: t_step,vmod,coll_freq
            logical, intent(in) :: boole_collisions,boole_random_precalc
            double precision, intent(in), optional :: nu_star
!
            integer :: i,j,k,l,n,o,nskip,ierr,iface,ind_tetr,iper
            integer :: counter_tetrahedron_passes, counter_lost_particles,counter_mappings
            integer :: n_lost_unit
            logical :: boole_t_finished
            logical, dimension(:), allocatable :: lost_particles
            integer, dimension(:), allocatable :: lost_reason
            integer(kind=8), dimension(:), allocatable :: lost_step
            integer, dimension(:), allocatable :: n_contributors
            double precision, dimension(:,:), allocatable :: lost_position
            double precision, dimension(:,:), allocatable :: psi_mat,psi_vec,delta_psi_mat,delta_psi_vec
            double precision, dimension(:,:), allocatable :: delta_psi2_average,delta_psi4_average,xi_collisions
            double precision, dimension(:,:), allocatable :: trend_delta_psi2_average,trend_std_delta_psi2_average
            double precision, dimension(:), allocatable :: xi
            double precision :: t_pass,t_remain,vpar,vperp,vmod1,pitchpar,collisionality,eps_collisions
            double precision, dimension(3) :: x_rand
!
            integer :: ipert
            double precision :: bmod_multiplier
            double precision :: A_s,A_theta,A_phi,B_s,B_theta,B_phi,bmod,q,dR_ds,dZ_ds,sqg
            double precision :: delta_A_amplitude,xi_collision
            integer,dimension(:), allocatable :: num_vec
            integer :: num_vec_size
            logical :: boole_initialized
            double precision, dimension(5) :: z
            double precision :: dtau, dphi,dtaumin,tau_out_can
            double precision :: theta_symflux
            double precision :: diff_coef, std_diff_coef, off_set
            double precision :: delta_psi
            integer(kind=8) :: fit_start_step,fit_end_step,valid_steps
            integer :: kpart = 0 ! progress counter for particles
            logical :: use_batch
            double precision, allocatable :: batch_x(:,:),batch_vpar(:),batch_vperp(:),batch_psi0(:)
            integer, allocatable :: batch_tetr(:),batch_iface(:)
            logical, allocatable :: batch_initialized(:)

            use_batch = trim(transport_backend) /= 'cpu'
            if(use_batch) then
                if(i_integrator_type /= 1) error stop 'batch transport requires i_integrator_type=1'
                if(.not.boole_random_precalc) error stop 'batch transport requires precomputed random numbers'
                if(boole_helical_pert) error stop 'batch transport does not support helical perturbations'
                allocate(batch_x(3,n_particles),batch_vpar(n_particles),batch_vperp(n_particles))
                allocate(batch_psi0(n_particles),batch_tetr(n_particles),batch_iface(n_particles))
                allocate(batch_initialized(n_particles))
            endif
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
                dtau = -1.d0*t_step*vmod
!
                !field line integration step step over phi (to check chamber wall crossing)
                dphi=2.d0*pi/dble(100)
!
                !orbit integration time step (to check chamber wall crossing)
                dtaumin=dphi*mag_axis_R0
!
                ierr = 0
            endif
!
!------------------------------------------------------------------------------------------------------------!
!---------------- Allocate and initialize variales for diffusion coefficient computation --------------------!
!
            if(boole_psi_mat) then
                allocate(psi_mat(n_particles,n_time_steps+1))
                psi_mat = 0.d0
                allocate(delta_psi_mat(n_particles,n_time_steps))
                delta_psi_mat = 0.d0
            else
                allocate(psi_vec(1,n_time_steps+1))
                allocate(delta_psi_vec(1,n_time_steps))
            endif
!
            allocate(delta_psi2_average(1,n_time_steps))
            allocate(delta_psi4_average(1,n_time_steps))
            allocate(trend_delta_psi2_average(2,n_time_steps))
            allocate(trend_std_delta_psi2_average(2,n_time_steps))
            allocate(xi(n_particles))
            allocate(lost_particles(n_particles))
            allocate(lost_reason(n_particles),lost_step(n_particles))
            allocate(lost_position(3,n_particles))
            allocate(n_contributors(n_time_steps))
!
            !Initializations
            delta_psi2_average = 0.d0
            delta_psi4_average = 0.d0
!
            !Particles are lost in MC Simulation
            lost_particles = .false. !Standard false, if particles are lost: true.
            lost_reason = 0
            lost_step = 0
            lost_position = 0.d0
            n_contributors = 0
            counter_lost_particles = 0
!
            !collisionality (Boozer & Kuo-Petravic 1981)
            eps_collisions = coll_freq*t_step
!
!------------------------------------------------------------------------------------------------------------!
!------------------------------------ Random numbers for MC simulation --------------------------------------!
!
            !Precompute random numbers
            if(boole_random_precalc) then
                call calc_rand_numbers(n_particles,n_time_steps)
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
!num_vec_size = 1
!allocate(num_vec(num_vec_size))
!num_vec = & ! [62,188,313,437,560,689,811,938,1062,1184,1313,1437] ! &
!       &  [2060] !,2188] !,2314,2439,2560,2691,2815,2936]
!

!
! open(89,file='error_data.dat')
! write(89,*) n_particles
!print *, n_particles
!
!------------------------------------------------------------------------------------------------------------!
!---------------------------------------- Loop over all particles -------------------------------------------!
!
            !$OMP PARALLEL DEFAULT(NONE) &
            !$OMP& SHARED(n_particles,pos_fluxtv_mat,xi,t_step,psi_mat,boole_random_precalc, &
            !$OMP& vmod,coord_system,n_time_steps,lost_particles,psitor_max,i_integrator_type,cm_over_e, &
            !$OMP& lost_reason,lost_step,lost_position, &
            !$OMP& n_contributors, &
            !$OMP& use_batch,batch_x,batch_vpar,batch_vperp,batch_psi0,batch_tetr,batch_iface,batch_initialized, &
            !$OMP& boole_collisions,eps_collisions,rd_start_pitchpar,rd_xi_collision,num_vec,kpart,boole_psi_mat, &
            !$OMP& delta_psi2_average, delta_psi4_average,lambda_start,helical_pert_eps_Aphi,helical_pert_m_fourier, &
            !$OMP& helical_pert_n_fourier,boole_helical_pert) &
            !$OMP& FIRSTPRIVATE(n,vpar,vperp,i,x_rand,ind_tetr,iface,k,z,psi_vec,delta_psi_vec, &
            !$OMP& counter_tetrahedron_passes,t_remain,dtau,dtaumin, &
            !$OMP& t_pass,boole_t_finished,iper,ierr, &
            !$OMP& ipert,bmod_multiplier,A_s,A_theta,A_phi,B_s,B_theta,B_phi,boole_initialized,dR_ds,dZ_ds,sqg, &
            !$OMP& bmod,q,delta_A_amplitude,counter_mappings,vmod1,pitchpar,xi_collision,l,num_vec_size,tau_out_can, &
            !$OMP& valid_steps)
            !$OMP DO
!
!do o = 1,num_vec_size
!n = num_vec(o)
!
            do n = 1,n_particles
!
                !Initialize psi_vec & delta_psi_vec
                if(.not.boole_psi_mat) then
                    psi_vec = 0.d0
                    delta_psi_vec = 0.d0
                endif
!
                !Counter for particles
                !$omp critical
                kpart = kpart+1
!print *, kpart, ' / ', n_particles, 'particle: ', n, 'thread: ', omp_get_thread_num()
                !$omp end critical
!
! print *, 'n',n
!
!if((i_integrator_type.eq.1).or.(i_integrator_type.eq.0)) then
!    !Particle velocities for every single Monte Carlo run
!    vpar = 1.d0*vmod
!    vperp = 0.d0!1.d-16*vmod
!endif
!
                !Pitchparameter at start - Random starting pitchangle in case of collisions
                if(boole_collisions) then
                    if(boole_random_precalc) then
                        pitchpar = rd_start_pitchpar(n)
                    else
                        !$omp critical
                            call random_number(pitchpar)
                        !$omp end critical
                    endif
                    pitchpar = 2.d0*pitchpar - 1.d0
!
                !Manually set pitchparameter for all particles
                else
                    pitchpar = lambda_start
                endif
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
                !Find random start indices that are distributed proportionally to the flux tube volume
                
                call binsrc(pos_fluxtv_mat(:,4),1,size(pos_fluxtv_mat(:,4)),xi(n),i)
!
                select case(i_integrator_type)
                    case(1,0)
                        x_rand = pos_fluxtv_mat(i,1:3)
                    case(2)
                        z(1) = pos_fluxtv_mat(i,1)
                        z(3) = pos_fluxtv_mat(i,3)
                        z(2) = theta_sym_flux2theta_vmec(z(1),pos_fluxtv_mat(i,2),z(3))  !Transform theta_symflux to theta_vmec
                        z(4) = 1.d0
                end select
!
!------------------------------------------------------------------------------------------------------------!
!--------------------------------- Option for an analytical perturbation ------------------------------------!
!
                !Manipulate psi starting value according to analytical perturbation function.
                !The s-label, which is not a flux surface label anymore in the perturbed system, is changed in linear approximation so that all particles
                !lie on the same flux surface of perturbed field
                if(boole_helical_pert.and.(coord_system.eq.2)) then
                    bmod_multiplier = 1.d0
                    ipert = 0
                    call vector_potential_sthetaphi(x_rand(1),x_rand(2),x_rand(3),ipert,bmod_multiplier,A_s,A_theta, &
                                                      & A_phi,B_s,B_theta,B_phi,bmod,q,dR_ds,dZ_ds)
!
                    !Manipulate s --> Compensate it so that s + delta_s is a flux label
                    delta_A_amplitude = helical_pert_m_fourier * A_phi*helical_pert_eps_Aphi / &
                                       & (helical_pert_n_fourier + helical_pert_m_fourier/q) / psitor_max
                endif
!
!------------------------------------------------------------------------------------------------------------!
!------------------- Initialization of integrators and initial psi measurement ------------------------------!
!
                select case(i_integrator_type)
                    case(0) !Modelling without orbit computation (speed offset)
                        !Initial Psi measurement
                        select case(coord_system)
                            case(1) !R,Phi,Z --> Cylindrical coordinate system
                                if(boole_psi_mat) then
                                    psi_mat(n,1) = psi_rphiz(ind_tetr,x_rand)
                                else
                                    psi_vec(1,1) = psi_rphiz(ind_tetr,x_rand)
                                endif
!
                            case(2) !s,theta,phi --> Symmetry flux coordinate system
                                !Save flux and subtract the bias due to the perturbation of the vector potential
                                if(boole_psi_mat) then
                                    psi_mat(n,1) = x_rand(1)
                                else
                                    psi_vec(1,1) = x_rand(1)
                                endif
                        end select
!
                    case(1)
!
                        !Initial Psi measurement
                        select case(coord_system)
                            case(1) !R,Phi,Z --> Cylindrical coordinate system
                               !Initialize particle in grid
                                boole_initialized = .false.
                                call orbit_timestep_gorilla(x_rand,vpar,vperp,0.d0,boole_initialized,ind_tetr,iface)

                                if(use_batch) then
                                    if(ind_tetr < 1) error stop 'batch transport start lies outside mesh'
                                endif

                                if(boole_psi_mat) then
                                    psi_mat(n,1) = psi_rphiz(ind_tetr,x_rand)
                                else
                                    psi_vec(1,1) = psi_rphiz(ind_tetr,x_rand)
                                endif
!
                            case(2) !s,theta,phi --> Symmetry flux coordinate system
                                !Save flux and subtract the bias due to the perturbation of the vector potential
                                if(boole_psi_mat) then
                                    ! Use corrugated flux surface label s_0
                                    psi_mat(n,1) = x_rand(1)
                                else
                                    ! Use corrugated flux surface label s_0
                                    psi_vec(1,1) = x_rand(1)
                                endif
!
                                if(boole_helical_pert) then
                                    ! Set particle manually to flux surface label s = s_0 + \delta s
                                    x_rand(1) = x_rand(1) + delta_A_amplitude  &
                                              & * cos(helical_pert_m_fourier*x_rand(2) +helical_pert_n_fourier*x_rand(3))
                                endif
!
                                !Initialize particle in grid
                                boole_initialized = .false.
                                call orbit_timestep_gorilla(x_rand,vpar,vperp,0.d0,boole_initialized,ind_tetr,iface)
                                if(use_batch) then
                                    if(ind_tetr < 1) error stop 'batch transport start lies outside mesh'
                                endif
                        end select
!
                    case(2)
                        !Initial psi measurement
                        if(boole_psi_mat) then
                            psi_mat(n,1) = z(1)
                        else
                            psi_vec(1,1) = z(1)
                        endif
                end select !i_integrator_type

                if(use_batch) then
                    batch_x(:,n) = x_rand
                    batch_vpar(n) = vpar
                    batch_vperp(n) = vperp
                    batch_tetr(n) = ind_tetr
                    batch_iface(n) = iface
                    batch_initialized(n) = boole_initialized
                    if(boole_psi_mat) then
                        batch_psi0(n) = psi_mat(n,1)
                    else
                        batch_psi0(n) = psi_vec(1,1)
                    endif
                    cycle
                endif
!
!------------------------------------------------------------------------------------------------------------!
!----------------------------------- Time steps of individual particles -------------------------------------!
!
                do k = 1,n_time_steps
!
! print *, 'n',n,'k',k
!if((n.eq.2060).and.(k.eq.152305)) transp_coef_trigger = .true.
!
                    !Orbit integration and measurement of psi after this time step
                    select case(i_integrator_type)
                        case(0)
                            !Psi measurement at x(t')
                            select case(coord_system)
                                case(1) !R,Phi,Z --> Cylindrical coordinate system
                                    if(boole_psi_mat) then
                                        psi_mat(n,k+1) = psi_rphiz(ind_tetr,x_rand)
                                    else
                                        psi_vec(1,k+1) = psi_rphiz(ind_tetr,x_rand)
                                    endif
!
                                case(2) !s,theta,phi --> Symmetry flux coordinate system
                                    if(boole_psi_mat) then
                                        psi_mat(n,k+1) = x_rand(1)
                                    else
                                        psi_vec(1,k+1) = x_rand(1)
                                    endif
                            end select
!
                        case(1)
                            call orbit_timestep_gorilla(x_rand,vpar,vperp,t_step,boole_initialized,ind_tetr,iface)
!
                            if(ind_tetr.eq.-1) then
                                !print *,"Error in tetra_main_orb. Particle has left the domain at ind_tetr:",ind_tetr
                                !print *,'Position of the particle',x_rand
                                lost_particles(n) = .true.
                                lost_reason(n) = 1
                                lost_step(n) = k
                                lost_position(:,n) = x_rand
                                exit
                            endif
!
                            !Psi measurement at x(t')
                            select case(coord_system)
                                case(1) !R,Phi,Z --> Cylindrical coordinate system
                                    if(boole_psi_mat) then
                                        psi_mat(n,k+1) = psi_rphiz(ind_tetr,x_rand)
                                    else
                                        psi_vec(1,k+1) = psi_rphiz(ind_tetr,x_rand)
                                    endif
!
                                case(2) !s,theta,phi --> Symmetry flux coordinate system
                                    if(boole_psi_mat) then
                                        if(boole_helical_pert) then
                                            ! Use corrugated flux surface label s_0
                                            psi_mat(n,k+1) = x_rand(1)  &
                                            & - delta_A_amplitude  &
                                            & * cos(helical_pert_m_fourier*x_rand(2) +helical_pert_n_fourier*x_rand(3))
                                        else
                                            psi_mat(n,k+1) = x_rand(1)
                                        endif
                                    else
                                        if(boole_helical_pert) then
                                            ! Use corrugated flux surface label s_0
                                            psi_vec(1,k+1) = x_rand(1)  &
                                            & - delta_A_amplitude  &
                                            & * cos(helical_pert_m_fourier*x_rand(2) +helical_pert_n_fourier*x_rand(3))
                                        else
                                            psi_vec(1,k+1) = x_rand(1)
                                        endif
                                    endif
                            end select
!
                        case(2)
                            call orbit_timestep_can(z,dtau,dtaumin,ierr,tau_out_can)
!
                            if(ierr.eq.1) then
                                lost_particles(n) = .true.
                                lost_reason(n) = 2
                                lost_step(n) = k
                                lost_position(:,n) = z(1:3)
                                exit
                            endif
!
                            !Psi measurement at x(t')
                            if(boole_psi_mat) then
                                psi_mat(n,k+1) = z(1)
                            else
                                psi_vec(1,k+1) = z(1)
                            endif
                    end select !
!
                    !Introduce collisions  (Boozer & Kuo-Petravic 1981)
                    if(boole_collisions) then
!
                        select case(i_integrator_type)
                            case(1,0)
                                !Calculation of pitchparameter
                                vmod1 = sqrt(vpar**2 + vperp**2)
                                pitchpar = vpar/vmod1
                            case(2)
                                pitchpar = z(5)
                        end select
!
                        if(boole_random_precalc) then
                            xi_collision = rd_xi_collision(n,k)
                        else
                            !$omp critical
                                call random_number(xi_collision)
                            !$omp end critical
                        endif
                        if(xi_collision.gt.0.5d0) then
                            xi_collision = 1.d0
                        else
                            xi_collision = -1.d0
                        endif
                        pitchpar = lorentz_pitch_step(pitchpar,eps_collisions,xi_collision)
!
                        select case(i_integrator_type)
                            case(1,0)
                                !Update the velocities after collision due to pitchangle scattering
                                vpar = vmod1*pitchpar
                                vperp = sqrt(vmod1**2 - vpar**2)
                            case(2)
                                z(5) = pitchpar
                        end select
                    endif
!
                enddo !k = 1,n_time_steps
!
!------------------------------------------------------------------------------------------------------------!
!----------------------------- Treatment in case of NOT storing psi mat ----------------------------------!
!
                !Calculation of momenta of Delta Psi (t)
                if(.not.boole_psi_mat) then
                    valid_steps = n_time_steps
                    if(lost_particles(n)) valid_steps = lost_step(n)-1
!
                    if(valid_steps > 0) then
                        !Compute Delta Psi through the last successful orbit step.
                        do i = 1,valid_steps
                            delta_psi_vec(1,i) = psi_vec(1,i+1) - psi_vec(1,1)
                        enddo
!
                        !Delta Psi squared
                        delta_psi_vec(1,1:valid_steps) = delta_psi_vec(1,1:valid_steps)**2
!
                        !Accumulate moments only through the last successful orbit step.
                        !$omp critical
                        call accumulate_valid_moments(delta_psi_vec(1,:),valid_steps, &
                                                      & delta_psi2_average(1,:),delta_psi4_average(1,:), &
                                                      & n_contributors)
                        !$omp end critical
                    endif

                endif
!
!------------------------------------------------------------------------------------------------------------!
!
            enddo !n_particles (Monte Carlo runs)
            !$OMP END DO
            !$OMP END PARALLEL

            if(use_batch) then
                call advance_transport_batch(batch_x,batch_vpar,batch_vperp,batch_initialized,batch_tetr, &
                    & batch_iface,batch_psi0,n_time_steps,t_step,boole_collisions,eps_collisions,psi_mat, &
                    & lost_particles,lost_reason,lost_step,lost_position,delta_psi2_average(1,:), &
                    & delta_psi4_average(1,:),n_contributors)
            endif
!
!------------------------------------------------------------------------------------------------------------!
!-------------------------- Calculation of momenta of Delta Psi (t) - ensemble average ----------------------!
!
            do n = 1,n_particles
                if(lost_particles(n)) counter_lost_particles = counter_lost_particles+1
            enddo

            !If psi mat is stored
            if(boole_psi_mat) then
                do n = 1,n_particles
                    valid_steps = n_time_steps
                    if(lost_particles(n)) valid_steps = lost_step(n)-1
                    do i = 1,valid_steps
                        delta_psi = psi_mat(n,i+1)-psi_mat(n,1)
                        delta_psi2_average(1,i) = delta_psi2_average(1,i)+delta_psi**2
                        delta_psi4_average(1,i) = delta_psi4_average(1,i)+delta_psi**4
                        n_contributors(i) = n_contributors(i)+1
                    enddo
                enddo
            endif
!
            !Compute moments with every trajectory retained through its last valid step.
            do i = 1,n_time_steps
                if(n_contributors(i) > 0) then
                    trend_delta_psi2_average(1,i) =  i*t_step
                    delta_psi2_average(1,i) = delta_psi2_average(1,i)/n_contributors(i)
                    trend_delta_psi2_average(2,i) =  delta_psi2_average(1,i)
    !
                    delta_psi4_average(1,i)  = delta_psi4_average(1,i)/n_contributors(i)
                    trend_std_delta_psi2_average(1,i) = i*t_step
                    trend_std_delta_psi2_average(2,i) = sqrt(max(0.d0,delta_psi4_average(1,i) &
                                                        & - delta_psi2_average(1,i)**2)) &
                                                        & /sqrt(dble(n_contributors(i)))
    !
                    if( (idiffcoef_output.eq.2).or.(idiffcoef_output.eq.3) ) then
                        write(file_id_psi2,*) trend_delta_psi2_average(1,i),trend_delta_psi2_average(2,i)
                        write(file_id_std_psi2,*) trend_std_delta_psi2_average(1,i),trend_std_delta_psi2_average(2,i)
                    endif
                endif
            enddo

            ! Write an zero line at the end of the calculation for one collisionality
            if( (idiffcoef_output.eq.2).or.(idiffcoef_output.eq.3) ) then
                write(file_id_psi2,*) 0,0
                write(file_id_std_psi2,*) 0,0
            endif

!
            print *, 'Number of lost particles in Monte Carlo simulation', counter_lost_particles
            if(present(file_id_loss_summary)) then
                if(present(nu_star)) then
                    write(file_id_loss_summary,*) nu_star,counter_lost_particles,n_particles
                else
                    write(file_id_loss_summary,*) counter_lost_particles,n_particles
                endif
            endif
            if(present(file_id_loss_events)) then
                do n = 1,n_particles
                    if(.not.lost_particles(n)) cycle
                    if(present(nu_star)) then
                        write(file_id_loss_events,*) nu_star,n,lost_reason(n),lost_step(n), &
                                                    & lost_step(n)*t_step,lost_position(:,n)
                    else
                        write(file_id_loss_events,*) n,lost_reason(n),lost_step(n), &
                                                    & lost_step(n)*t_step,lost_position(:,n)
                    endif
                enddo
            endif
            if(present(file_id_particle_histories)) then
                if(.not.boole_psi_mat) then
                    error stop 'particle history output requires boole_psi_mat'
                endif
                do n = 1,n_particles
                    valid_steps = n_time_steps
                    if(lost_particles(n)) valid_steps = lost_step(n)-1
                    do i = 0,valid_steps
                        if(present(nu_star)) then
                            write(file_id_particle_histories,*) nu_star,n,i,i*t_step,psi_mat(n,i+1)
                        else
                            write(file_id_particle_histories,*) n,i,i*t_step,psi_mat(n,i+1)
                        endif
                    enddo
                enddo
            endif
!
            open(newunit=n_lost_unit,file='n_lost_particles.dat',status='unknown')
            write(n_lost_unit,*) counter_lost_particles
            close(n_lost_unit)
!
!------------------------------------------------------------------------------------------------------------!
!-------------- Computation of transport diffusion coefficient by linear least squares fit ------------------!
!
            if( (idiffcoef_output.eq.1).or.(idiffcoef_output.eq.3) ) then
!
                fit_start_step = ceiling(dble(n_time_steps)*0.2d0)
                fit_end_step = n_time_steps
                do n = 1,n_particles
                    if(.not.lost_particles(n)) cycle
                    fit_end_step = min(fit_end_step,lost_step(n)-1)
                enddo
                if(fit_end_step >= fit_start_step+1) then
                    if(fit_end_step < n_time_steps) then
                        print *, 'MSD fit stops before first particle loss at step',fit_end_step+1
                    endif
!
                    call llsq (fit_end_step-fit_start_step+1, &
                               & trend_delta_psi2_average(1,fit_start_step:fit_end_step), &
                               & trend_delta_psi2_average(2,fit_start_step:fit_end_step),diff_coef,off_set)
                    diff_coef = diff_coef/2.d0
!
                    call llsq (fit_end_step-fit_start_step+1, &
                               & trend_std_delta_psi2_average(1,fit_start_step:fit_end_step), &
                               & trend_std_delta_psi2_average(2,fit_start_step:fit_end_step),std_diff_coef,off_set)
                    std_diff_coef = std_diff_coef/2.d0
                else
                    print *, 'WARNING: no loss-free MSD fitting interval after transient, writing NaN'
                    diff_coef = ieee_value(diff_coef,ieee_quiet_nan)
                    std_diff_coef = ieee_value(std_diff_coef,ieee_quiet_nan)
                endif
!
                !Write output depending, if nu_star is present
                if(present(nu_star)) then
                    write(file_id_transp_diff_coef,*) nu_star, diff_coef, std_diff_coef
                else
                    write(file_id_transp_diff_coef,*) diff_coef, std_diff_coef
                endif
!
            endif
!
!------------------------------------------------------------------------------------------------------------!
!
            deallocate(delta_psi2_average,delta_psi4_average,xi,lost_particles)
            deallocate(lost_reason,lost_step,lost_position)
            deallocate(n_contributors)
            deallocate(trend_delta_psi2_average,trend_std_delta_psi2_average)
            if(boole_psi_mat) then
                deallocate(psi_mat,delta_psi_mat)
            else
                deallocate(psi_vec,delta_psi_vec)
            endif
!
            !Deallocate random numbers
            if(boole_random_precalc) then
                call dealloc_rand_numbers()
            endif
!
        end subroutine calc_flux_deviation

        subroutine advance_transport_batch(x,vpar,vperp,initialized,ind_tetr,iface,psi0,n_time_steps,t_step, &
                & boole_collisions,eps_collisions,psi_mat,lost,lost_reason,lost_step,lost_position, &
                & sum_psi2,sum_psi4,n_contributors)

            use orbit_timestep_gorilla_batch_mod, only: gorilla_batch_t,gorilla_batch_profile_t,initialize_gorilla_batch, &
                & finalize_gorilla_batch,orbit_timestep_gorilla_batch
            use gorilla_status_mod, only: GORILLA_COMPLETE,GORILLA_PHYSICAL_LOSS
            use mono_energetic_transp_coef_settings_mod, only: transport_backend
            use tetra_physics_mod, only: coord_system,psi_rphiz
            use omp_lib, only: omp_get_wtime
            use, intrinsic :: iso_fortran_env, only: int64

            implicit none

            double precision, intent(inout) :: x(:,:),vpar(:),vperp(:)
            logical, intent(inout) :: initialized(:)
            integer, intent(inout) :: ind_tetr(:),iface(:)
            double precision, intent(in) :: psi0(:),t_step,eps_collisions
            integer(kind=8), intent(in) :: n_time_steps
            logical, intent(in) :: boole_collisions
            double precision, allocatable, intent(inout) :: psi_mat(:,:)
            logical, intent(inout) :: lost(:)
            integer, intent(inout) :: lost_reason(:),n_contributors(:)
            integer(kind=8), intent(inout) :: lost_step(:)
            double precision, intent(inout) :: lost_position(:,:),sum_psi2(:),sum_psi4(:)

            type(gorilla_batch_t) :: batch
            type(gorilla_batch_profile_t) :: profile
            double precision, allocatable :: active_x(:,:),active_vpar(:),active_vperp(:)
            integer, allocatable :: active_ids(:),active_tetr(:),active_iface(:),status(:)
            logical, allocatable :: active_initialized(:)
            integer :: n,n_active,j,ierr,n_particles
            integer(kind=8) :: k
            integer(int64) :: orbit_advances,cell_crossings,replay_cells
            double precision :: psi,delta_psi,vmod,pitch,xi_collision,wall_start
            double precision :: setup_seconds,upload_seconds,kernel_seconds,download_seconds,replay_seconds
            character(16) :: backend

            n_particles = size(vpar)
            allocate(active_x(3,n_particles),active_vpar(n_particles),active_vperp(n_particles))
            allocate(active_ids(n_particles),active_tetr(n_particles),active_iface(n_particles))
            allocate(active_initialized(n_particles),status(n_particles))
            backend = transport_backend
            if(backend == 'cpu-batch') backend = 'cpu'
            wall_start = omp_get_wtime()
            call initialize_gorilla_batch(batch,trim(backend),ierr)
            if(ierr /= 0) then
                call finalize_gorilla_batch(batch)
                print *, 'Cannot initialize transport backend ',trim(transport_backend),'; error ',ierr
                error stop 'transport batch initialization failed'
            endif

            orbit_advances = 0_int64
            cell_crossings = 0_int64
            replay_cells = 0_int64
            setup_seconds = 0.d0
            upload_seconds = 0.d0
            kernel_seconds = 0.d0
            download_seconds = 0.d0
            replay_seconds = 0.d0
            do k = 1,n_time_steps
                n_active = 0
                do n = 1,n_particles
                    if(lost(n)) cycle
                    n_active = n_active+1
                    active_ids(n_active) = n
                    active_x(:,n_active) = x(:,n)
                    active_vpar(n_active) = vpar(n)
                    active_vperp(n_active) = vperp(n)
                    active_initialized(n_active) = initialized(n)
                    active_tetr(n_active) = ind_tetr(n)
                    active_iface(n_active) = iface(n)
                enddo
                if(n_active == 0) exit
                orbit_advances = orbit_advances+int(n_active,int64)
                call orbit_timestep_gorilla_batch(batch,active_x(:,1:n_active),active_vpar(1:n_active), &
                    & active_vperp(1:n_active),t_step,active_initialized(1:n_active),active_tetr(1:n_active), &
                    & active_iface(1:n_active),status(1:n_active),profile,ierr)
                if(ierr /= 0) then
                    call finalize_gorilla_batch(batch)
                    print *, 'Transport batch error ',ierr,', step ',k
                    error stop 'transport batch advance failed'
                endif
                cell_crossings = cell_crossings+profile%cell_crossings
                replay_cells = replay_cells+profile%replay_cells
                setup_seconds = setup_seconds+profile%setup_seconds
                upload_seconds = upload_seconds+profile%upload_seconds
                kernel_seconds = kernel_seconds+profile%kernel_seconds
                download_seconds = download_seconds+profile%download_seconds
                replay_seconds = replay_seconds+profile%replay_seconds

                do j = 1,n_active
                    n = active_ids(j)
                    x(:,n) = active_x(:,j)
                    vpar(n) = active_vpar(j)
                    vperp(n) = active_vperp(j)
                    initialized(n) = active_initialized(j)
                    ind_tetr(n) = active_tetr(j)
                    iface(n) = active_iface(j)
                    select case(status(j))
                        case(GORILLA_PHYSICAL_LOSS)
                            lost(n) = .true.
                            lost_reason(n) = 1
                            lost_step(n) = k
                            lost_position(:,n) = x(:,n)
                            cycle
                        case(GORILLA_COMPLETE)
                            if(ind_tetr(n) < 1) then
                                call finalize_gorilla_batch(batch)
                                error stop 'completed transport marker has no tetrahedron'
                            endif
                        case default
                            call finalize_gorilla_batch(batch)
                            print *, 'Transport numerical failure: marker ',n,', step ',k,', status ',status(j)
                            error stop 'numerical failure prevents a complete transport estimate'
                    end select

                    select case(coord_system)
                        case(1)
                            psi = psi_rphiz(ind_tetr(n),x(:,n))
                        case(2)
                            psi = x(1,n)
                        case default
                            call finalize_gorilla_batch(batch)
                            error stop 'unsupported transport coordinate system'
                    end select
                    if(allocated(psi_mat)) then
                        psi_mat(n,k+1) = psi
                    else
                        delta_psi = psi-psi0(n)
                        call accumulate_valid_moments([delta_psi**2],1_8,sum_psi2(k:k), &
                            & sum_psi4(k:k),n_contributors(k:k))
                    endif

                    if(boole_collisions) then
                        vmod = sqrt(vpar(n)**2+vperp(n)**2)
                        pitch = vpar(n)/vmod
                        xi_collision = -1.d0
                        if(rd_xi_collision(n,k) > 0.5d0) xi_collision = 1.d0
                        pitch = lorentz_pitch_step(pitch,eps_collisions,xi_collision)
                        vpar(n) = vmod*pitch
                        vperp(n) = sqrt(vmod**2-vpar(n)**2)
                    endif
                enddo
            enddo
            call finalize_gorilla_batch(batch)
            print *, 'Transport batch backend ',trim(transport_backend)
            print *, 'Transport batch markers ',n_particles,', completed ',count(.not.lost),', physical losses ',count(lost)
            print *, 'Transport batch numerical failures ',0
            print *, 'Transport batch orbit advances ',orbit_advances
            print *, 'Transport batch cell crossings ',cell_crossings,', replay cells ',replay_cells
            print *, 'Transport batch setup time [s] ',setup_seconds
            print *, 'Transport batch upload time [s] ',upload_seconds
            print *, 'Transport batch kernel time [s] ',kernel_seconds
            print *, 'Transport batch download time [s] ',download_seconds
            print *, 'Transport batch replay time [s] ',replay_seconds
            print *, 'Transport batch wall time [s] ',omp_get_wtime()-wall_start

        end subroutine advance_transport_batch
!
    end module flux_deviation_mod
