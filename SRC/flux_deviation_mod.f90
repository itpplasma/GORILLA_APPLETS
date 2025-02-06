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
        public :: calc_flux_deviation
!
    contains
!
        subroutine calc_rand_numbers(n_particles,n_time_steps)
!
            use mono_energetic_transp_coef_settings_mod, only: seed_option
!
            implicit none
!
            integer, intent(in)     :: n_particles
            integer(kind=8),intent(in) :: n_time_steps
            integer,dimension(:), allocatable :: seed
            double precision, dimension(:), allocatable :: rd_seed
            integer :: j,n
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
                            & nu_star)
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
            use mono_energetic_transp_coef_settings_mod, only: lambda_start,boole_psi_mat
            use gorilla_settings_mod, only: boole_helical_pert, helical_pert_eps_Aphi, helical_pert_m_fourier, &
                                          & helical_pert_n_fourier
            use binsrc_mod, only: binsrc
!
            implicit none
!
            integer,intent(in) :: n_particles,file_id_psi2,file_id_std_psi2,i_integrator_type,file_id_transp_diff_coef
            integer,intent(in) :: idiffcoef_output
            integer(kind=8),intent(in) :: n_time_steps
            double precision, intent(in) :: t_step,vmod,coll_freq
            logical, intent(in) :: boole_collisions,boole_random_precalc
            double precision, intent(in), optional :: nu_star
!
            integer :: i,j,k,l,n,o,nskip,ierr,iface,ind_tetr,iper
            integer :: counter_tetrahedron_passes, counter_lost_particles,counter_mappings
            logical :: boole_t_finished
            logical, dimension(:), allocatable :: lost_particles
            double precision, dimension(:,:), allocatable :: psi_mat,psi_vec,delta_psi_mat,delta_psi_vec
            double precision, dimension(:,:), allocatable :: delta_psi_average,delta_psi2_average,delta_psi4_average,xi_collisions
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
            integer :: kpart = 0 ! progress counter for particles
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
            allocate(delta_psi_average(1,n_time_steps))
            allocate(delta_psi2_average(1,n_time_steps))
            allocate(delta_psi4_average(1,n_time_steps))
            allocate(trend_delta_psi2_average(2,n_time_steps))
            allocate(trend_std_delta_psi2_average(2,n_time_steps))
            allocate(xi(n_particles))
            allocate(lost_particles(n_particles))
!
            !Initializations
            delta_psi_average = 0.d0
            delta_psi2_average = 0.d0
            delta_psi4_average = 0.d0
!
            !Particles are lost in MC Simulation
            lost_particles = .false. !Standard false, if particles are lost: true.
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
            !$OMP& boole_collisions,eps_collisions,rd_start_pitchpar,rd_xi_collision,num_vec,kpart,boole_psi_mat, &
            !$OMP& delta_psi2_average, delta_psi4_average,lambda_start,helical_pert_eps_Aphi,helical_pert_m_fourier, &
            !$OMP& helical_pert_n_fourier,boole_helical_pert) &
            !$OMP& FIRSTPRIVATE(n,vpar,vperp,i,x_rand,ind_tetr,iface,k,z,psi_vec,delta_psi_vec, &
            !$OMP& counter_tetrahedron_passes,t_remain,dtau,dtaumin, &
            !$OMP& t_pass,boole_t_finished,iper,ierr, &
            !$OMP& ipert,bmod_multiplier,A_s,A_theta,A_phi,B_s,B_theta,B_phi,boole_initialized,dR_ds,dZ_ds,sqg, &
            !$OMP& bmod,q,delta_A_amplitude,counter_mappings,vmod1,pitchpar,xi_collision,l,num_vec_size,tau_out_can)
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
                        pitchpar = pitchpar * (1.d0 - eps_collisions) + xi_collision * sqrt((1.d0-pitchpar**2)*eps_collisions)
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
!----------------------------- Treatment in case of of NOT storing psi mat ----------------------------------!
!
                !Calculation of momenta of Delta Psi (t)
                if(.not.boole_psi_mat) then
                    if(lost_particles(n)) cycle !If particles is lost, don't add up to the moments
!
                    !Compute Delta Psi
                    do i = 1,n_time_steps
                        delta_psi_vec(1,i) = psi_vec(1,i+1) - psi_vec(1,1)
                    enddo
!
                    !Delta Psi squared
                    delta_psi_vec = delta_psi_vec**2 !attention: Same variable is used
!
                    !Average Delta Psi Squared over all particles (only summation at this point)
                    !$omp critical
                    delta_psi2_average(1,:) = delta_psi2_average(1,:) + delta_psi_vec(1,:)
                    !$omp end critical
!
                    !Delta Psi to the power of four
                    delta_psi_vec = delta_psi_vec**2 !attention: Same variable is used

                    !Average Delta Psi to the power of four over all particles (only summation at this point)
                    !$omp critical
                    delta_psi4_average(1,:) = delta_psi4_average(1,:) + delta_psi_vec(1,:)
                    !$omp end critical

                endif
!
!------------------------------------------------------------------------------------------------------------!
!
            enddo !n_particles (Monte Carlo runs)
            !$OMP END DO
            !$OMP END PARALLEL
!
!------------------------------------------------------------------------------------------------------------!
!-------------------------- Calculation of momenta of Delta Psi (t) - ensemble average ----------------------!
!
            !If psi mat is stored
            if(boole_psi_mat) then
                !Delta Psi as a function of time for every single particle
                do i = 1,n_time_steps
                    delta_psi_mat(:,i) = psi_mat(:,i+1) - psi_mat(:,1)
                enddo
!
                !Exclude those particles, that have left the domain
                do i = 1,n_particles
                if(lost_particles(i)) then
                    delta_psi_mat(i,:) = 0.d0
                    counter_lost_particles = counter_lost_particles +1
                endif
                enddo
!
                !Delta Psi averaged over all particles
                delta_psi_average(1,:) = sum(delta_psi_mat,DIM=1)
                do i = 1,n_time_steps
                delta_psi_average(1,i)  = delta_psi_average(1,i)/(n_particles-counter_lost_particles)
                enddo
!
                !Delta Psi Squared averaged over all particles
                delta_psi_mat = delta_psi_mat**2 !attention: Same variable is used
                delta_psi2_average(1,:) = sum(delta_psi_mat,DIM=1)
!
                !Standard deviation of Delta Psi Squared
                delta_psi_mat = delta_psi_mat**2
                delta_psi4_average(1,:) = sum(delta_psi_mat,DIM=1)
!
            !If psi values are only stored for a single particle
            else
!
                !Exclude those particles, that have left the domain
                do i = 1,n_particles
                    if(lost_particles(i)) then
                        counter_lost_particles = counter_lost_particles +1
                    endif
                enddo
!
            endif
!
            !Compute momenta
            do i = 1,n_time_steps
                trend_delta_psi2_average(1,i) =  i*t_step
                delta_psi2_average(1,i) = delta_psi2_average(1,i)/(n_particles-counter_lost_particles)
                trend_delta_psi2_average(2,i) =  delta_psi2_average(1,i)
!
                delta_psi4_average(1,i)  = delta_psi4_average(1,i)/(n_particles-counter_lost_particles)
                trend_std_delta_psi2_average(1,i) = i*t_step
                trend_std_delta_psi2_average(2,i) = sqrt(delta_psi4_average(1,i) &
                                                    & - delta_psi2_average(1,i)**2)/sqrt(dble(n_particles))
!
                if( (idiffcoef_output.eq.2).or.(idiffcoef_output.eq.3) ) then
                    write(file_id_psi2,*) trend_delta_psi2_average(1,i),trend_delta_psi2_average(2,i)
                    write(file_id_std_psi2,*) trend_std_delta_psi2_average(1,i),trend_std_delta_psi2_average(2,i)
                endif
            enddo

            ! Write an zero line at the end of the calculation for one collisionality
            if( (idiffcoef_output.eq.2).or.(idiffcoef_output.eq.3) ) then
                write(file_id_psi2,*) 0,0
                write(file_id_std_psi2,*) 0,0
            endif

!
            print *, 'Number of lost particles in Monte Carlo simulation', counter_lost_particles
!
            open(unit=202,file='n_lost_particles.dat',status='unknown')
            write(202,*) counter_lost_particles
            close(202)
!
!------------------------------------------------------------------------------------------------------------!
!-------------- Computation of transport diffusion coefficient by linear least squares fit ------------------!
!
            if( (idiffcoef_output.eq.1).or.(idiffcoef_output.eq.3) ) then
!
                ! Starting index (Throw away first 20 percent of values)
                i = ceiling(dble(n_time_steps)*0.2d0)
!
                call llsq (n_time_steps-i+1 , trend_delta_psi2_average(1,i:n_time_steps) , &
                           &  trend_delta_psi2_average(2,i:n_time_steps) ,diff_coef,off_set )
                diff_coef = diff_coef/2.d0
!
                call llsq (n_time_steps-i+1 , trend_std_delta_psi2_average(1,i:n_time_steps) , &
                           &  trend_std_delta_psi2_average(2,i:n_time_steps) ,std_diff_coef,off_set )
                std_diff_coef = std_diff_coef/2.d0
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
            deallocate(delta_psi_average,delta_psi2_average,delta_psi4_average,xi,lost_particles)
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
!
    end module flux_deviation_mod
