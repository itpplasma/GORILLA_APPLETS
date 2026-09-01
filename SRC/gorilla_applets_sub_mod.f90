!
    module gorilla_applets_sub_mod
!
        implicit none
!
        private
!
        double precision            :: vmod
!
        public                      :: initialize_mono_energetic_transp_coef,calc_mono_energetic_transp_coef, &
                                    & calc_numerical_diff_coef,calc_mono_energetic_transp_coef_nu_scan, &
                                    & collision_frequency_from_nu_star,nu_star_from_collision_frequency, &
                                    & radial_metric_ds_dr
!
    contains
!
        pure function collision_frequency_from_nu_star(nu_star,major_radius,aiota,speed) result(nu_collision)
!
            implicit none
!
            double precision, intent(in) :: nu_star,major_radius,aiota,speed
            double precision :: nu_collision
!
            nu_collision = nu_star*abs(aiota)*speed/major_radius
!
        end function collision_frequency_from_nu_star
!
        pure function nu_star_from_collision_frequency(nu_collision,major_radius,aiota,speed) result(nu_star)
!
            implicit none
!
            double precision, intent(in) :: nu_collision,major_radius,aiota,speed
            double precision :: nu_star
!
            nu_star = major_radius*nu_collision/(abs(aiota)*speed)
!
        end function nu_star_from_collision_frequency
!
        pure function radial_metric_ds_dr(s,bmod,torflux_signed) result(ds_dr)
!
            use constants, only: pi
!
            implicit none
!
            double precision, intent(in) :: s,bmod,torflux_signed
            double precision :: ds_dr
!
            ds_dr = 2.d0*sqrt(s)*sqrt(pi*bmod/abs(torflux_signed))
!
        end function radial_metric_ds_dr
!
        subroutine initialize_mono_energetic_transp_coef()
!
            use orbit_timestep_gorilla_mod, only: initialize_gorilla
            use gorilla_settings_mod, only: set_eps_Phi
            use constants, only: clight,ev2erg,pi
            use tetra_physics_mod, only: particle_mass
            use tetra_grid_settings_mod, only: grid_kind
            use vector_potential_mod, only: torflux
            use mono_energetic_transp_coef_settings_mod, only: energy_eV,v_E
            use fluxtv_mod, only: pos_fluxtv_mat,load_flux_tube_volume
            use splint_vmec_data_mod, only: splint_vmec_data
            use sub_alpha_lifetime_can_mod, only: integrate_mfl_can
!
            implicit none
!
            double precision, dimension(3)              :: x_start
            double precision                            :: s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota
            double precision                            :: R,R0,Z,Z0,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
            double precision                            :: sqg,Bctrvr_vartheta,Bctrvr_varphi
            double precision                            :: Bcovar_r,Bcovar_vartheta,Bcovar_varphi,bmod,ds_dr
            double precision                            :: xi,mag_axis_bmod
            double precision                            :: eps_Phi
            integer                                     :: i,ierr
            integer,parameter :: npoi = 100000
            double precision :: dphi,phibeg,rbeg,zbeg,bmod00
            double precision, dimension(3,npoi) :: xstart
            double precision, dimension(npoi)   :: bstart,volstart
!
            call initialize_gorilla(1)
!
            !Compute velocity module from kinetic energy dependent on particle species
            vmod=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
            print *, 'vmod', vmod
!
            !Load fluxtube volume for a starting position (File to be chosen in gorilla_applets.inp)
            call load_flux_tube_volume()
            x_start = pos_fluxtv_mat(1,1:3)
!
            !----------- Compute electrostatic potential for normalized ExB drift velocity v_E -----------!
            if(grid_kind.eq.3) then
!
               dphi = 2.d0*pi/100.d0
               rbeg = x_start(1)
               phibeg = x_start(3)
               zbeg = x_start(2)
!
               print *, 'Start to integrate MFL to obtain bmod00.'
               call integrate_mfl_can(npoi,dphi,rbeg,phibeg,zbeg,         &
                                   xstart,bstart,volstart,bmod00,ierr)
               print *, 'Finished to integrate MFL to obtain bmod00.'
!
               !Compute E_r = dPhi_dr = dPhi_ds * ds_dr = dA_theta_ds * ds_dr
               s = x_start(1)
               theta = x_start(2)
               varphi = x_start(3)
!
               call splint_vmec_data(s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                                     R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
               print *, 'bmod',bmod00
!
               !Explanation for ds_dr: r = sqrt(s) * sqrt(Psi_tor_a/(pi * B00))
               ds_dr = radial_metric_ds_dr(s,bmod00,torflux)
!
               !Set electrostatic potential to values such that normalized elecron drift velocity is reproduced
               eps_Phi = v_E*vmod*bmod00/(dA_theta_ds*ds_dr*clight)!
               call set_eps_Phi(eps_Phi)
               print *, 'eps_Phi', eps_Phi!
!
           endif
!
           call initialize_gorilla(2)
!
        end subroutine initialize_mono_energetic_transp_coef
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine calc_mono_energetic_transp_coef()
!
            use constants, only: pi
            use tetra_physics_mod, only: mag_axis_R0,tetra_physics
            use tetra_grid_mod, only: ntetr
            use tetra_grid_settings_mod, only: grid_kind,grid_size,n_field_periods
            use flux_deviation_mod, only: calc_flux_deviation
            use mono_energetic_transp_coef_settings_mod, only: boole_collisions, boole_random_precalc, i_integrator_type, &
                                       & filename_transp_diff_coef, filename_delta_s_squared, filename_std_dvt_delta_s_squared, &
                                       & idiffcoef_output, energy_eV, v_E, nu_star, n_particles,flight_time_multiplier, &
                                       & boole_write_particle_histories,filename_particle_histories, &
                                       & filename_transport_metadata,seed_option,random_seed_filename
!
            implicit none
!
            integer                                     :: file_id_psi2,file_id_std_psi2,file_id_transp_diff_coef
            integer                                     :: file_id_loss_events,file_id_loss_summary
            integer                                     :: file_id_particle_histories,file_id_transport_metadata
            integer(kind=8)                             :: n_time_steps
            double precision                            :: tau_bounce,coll_freq,tau_collision,t_step
            double precision                            :: q_saf,aiota
!
            !Initialize GORILLA with electrostatic potential in accordance with Mach number for mono-energetic transport coefficient
            call initialize_mono_energetic_transp_coef()
!
            call get_local_rotational_transform(q_saf,aiota)
            print *, 'Local safety factor q', q_saf
            print *, 'Local rotational transform iota', aiota
!
            select case(idiffcoef_output)
                case(1)
                    open(newunit=file_id_transp_diff_coef,file=filename_transp_diff_coef)
                case(2)
                    open(newunit=file_id_psi2,file=filename_delta_s_squared)
                    open(newunit=file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
                case(3)
                    open(newunit=file_id_transp_diff_coef,file=filename_transp_diff_coef)
                    open(newunit=file_id_psi2,file=filename_delta_s_squared)
                    open(newunit=file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
            end select
            open(newunit=file_id_loss_events,file='lost_particle_events.dat',status='replace')
            open(newunit=file_id_loss_summary,file='loss_summary.dat',status='replace')
            open(newunit=file_id_transport_metadata,file=trim(filename_transport_metadata),status='replace')
            call write_transport_metadata_header(file_id_transport_metadata,random_seed_filename)
            if(boole_write_particle_histories) then
                open(newunit=file_id_particle_histories,file=trim(filename_particle_histories),status='replace')
            endif
!
            print *, 'Mono-energetic radial transport coefficient:'
            print *, ''
            print *, 'Normalized collisionality', nu_star
            print *, 'Mach number', v_E
            print *, 'Number of particles', n_particles
!
            select case(grid_kind)
                case(1,2,6)
                    !Define bounce time for Tokamak
                    !With the local q this is nu_star = R_0 nu_c q / v_mod = R_0 nu_c / (iota v_mod)
                    tau_bounce = 2.d0*pi*mag_axis_R0/vmod*abs(q_saf)
                    tau_collision = 1.d0/collision_frequency_from_nu_star(nu_star,mag_axis_R0,aiota,vmod)
                case(3)
                    !Define bounce time for Stellarator
                    tau_bounce = 2.d0*pi*mag_axis_R0/vmod/n_field_periods
                    tau_collision = 1.d0/collision_frequency_from_nu_star(nu_star,mag_axis_R0,aiota,vmod)
            end select
!
            !Define collision frequency from collision time
            coll_freq = 1.d0/tau_collision
!
            !Define step size
            t_step = minval([tau_bounce/20.d0,tau_collision/20.d0]);
!
            !Define number of steps
            n_time_steps = int(flight_time_multiplier &
                               & * maxval([tau_collision,tau_bounce**2/tau_collision])/t_step)

            call write_transport_metadata_row(file_id_transport_metadata,nu_star,mag_axis_R0,aiota,q_saf, &
                                               & vmod,coll_freq,t_step,n_time_steps,n_particles,seed_option, &
                                               & n_field_periods)

            print *, 'Number of time steps (measurements)',n_time_steps
!
            print *, 'Start calculation of diffusion coefficient'
!
            if(boole_write_particle_histories) then
                call calc_flux_deviation(n_particles,n_time_steps,t_step,vmod,boole_collisions,boole_random_precalc,coll_freq, &
                                    & file_id_transp_diff_coef,file_id_psi2,file_id_std_psi2,i_integrator_type,idiffcoef_output, &
                                    & nu_star,file_id_loss_events,file_id_loss_summary,file_id_particle_histories)
                close(file_id_particle_histories)
            else
                call calc_flux_deviation(n_particles,n_time_steps,t_step,vmod,boole_collisions,boole_random_precalc,coll_freq, &
                                    & file_id_transp_diff_coef,file_id_psi2,file_id_std_psi2,i_integrator_type,idiffcoef_output, &
                                    & nu_star,file_id_loss_events,file_id_loss_summary)
            endif
            close(file_id_loss_events)
            close(file_id_loss_summary)
            close(file_id_transport_metadata)
!
            select case(idiffcoef_output)
                case(1)
                    close(file_id_transp_diff_coef)
                case(2)
                    close(file_id_psi2)
                    close(file_id_std_psi2)
                case(3)
                    close(file_id_transp_diff_coef)
                    close(file_id_psi2)
                    close(file_id_std_psi2)
            end select
!
        end subroutine calc_mono_energetic_transp_coef
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine calc_mono_energetic_transp_coef_nu_scan()
!
            use constants, only: pi
            use tetra_physics_mod, only: mag_axis_R0,tetra_physics
            use tetra_grid_mod, only: ntetr
            use tetra_grid_settings_mod, only: grid_kind,grid_size,n_field_periods
            use flux_deviation_mod, only: calc_flux_deviation
            use mono_energetic_transp_coef_settings_mod, only: boole_collisions, boole_random_precalc, i_integrator_type, &
                                       & filename_transp_diff_coef, filename_delta_s_squared, filename_std_dvt_delta_s_squared, &
                                       & idiffcoef_output, energy_eV, v_E, n_particles, n_nu_scans, nu_star_start, &
                                       & nu_exp_basis,flight_time_multiplier,boole_write_particle_histories, &
                                       & filename_particle_histories,filename_transport_metadata,seed_option, &
                                       & random_seed_filename
!
            implicit none
!
            integer                                     :: i,file_id_psi2,file_id_std_psi2,file_id_transp_diff_coef
            integer                                     :: file_id_loss_events,file_id_loss_summary
            integer                                     :: file_id_particle_histories,file_id_transport_metadata
            integer(kind=8)                             :: n_time_steps
            double precision                            :: tau_bounce,coll_freq,tau_collision,t_step, nu_star
            double precision                            :: q_saf,aiota
!
            !Initialize GORILLA with electrostatic potential in accordance with Mach number for mono-energetic transport coefficient
            call initialize_mono_energetic_transp_coef()
!
            !q and iota are flux functions, so one evaluation ahead of the nu* loop is enough.
            call get_local_rotational_transform(q_saf,aiota)
            print *, 'Local safety factor q', q_saf
            print *, 'Local rotational transform iota', aiota
!
            select case(idiffcoef_output)
                case(1)
                    open(newunit=file_id_transp_diff_coef,file=filename_transp_diff_coef)
                case(2)
                    open(newunit=file_id_psi2,file=filename_delta_s_squared)
                    open(newunit=file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
                case(3)
                    open(newunit=file_id_transp_diff_coef,file=filename_transp_diff_coef)
                    open(newunit=file_id_psi2,file=filename_delta_s_squared)
                    open(newunit=file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
            end select
            open(newunit=file_id_loss_events,file='lost_particle_events.dat',status='replace')
            open(newunit=file_id_loss_summary,file='loss_summary.dat',status='replace')
            open(newunit=file_id_transport_metadata,file=trim(filename_transport_metadata),status='replace')
            call write_transport_metadata_header(file_id_transport_metadata,random_seed_filename)
            if(boole_write_particle_histories) then
                open(newunit=file_id_particle_histories,file=trim(filename_particle_histories),status='replace')
            endif
!
            print *, 'Mono-energetic radial transport coefficient - Scan over normalized collisionality:'
            print *, ''
            print *, 'Mach number', v_E
            print *, 'Number of particles', n_particles
            print *, ''
!
            !Scan over collisionalities
            do i = 0,n_nu_scans-1
!
                !Definition of nu_star
                nu_star = nu_star_start * nu_exp_basis**i
!
                print *, ''
                print *, 'Scan over normalized collisionality - Nr. ',i+1
                print *, 'Normalized collisionality', nu_star
                select case(grid_kind)
                    case(1,2,6)
                        !Define bounce time for Tokamak
                        tau_bounce = 2.d0*pi*mag_axis_R0/vmod*abs(q_saf)
                        tau_collision = 1.d0/collision_frequency_from_nu_star(nu_star,mag_axis_R0,aiota,vmod)
                    case(3)
                        !Define bounce time for Stellarator
                        tau_bounce = 2.d0*pi*mag_axis_R0/vmod/n_field_periods
                        tau_collision = 1.d0/collision_frequency_from_nu_star(nu_star,mag_axis_R0,aiota,vmod)
                end select
!
                !Define collision frequency from collision time
                coll_freq = 1.d0/tau_collision
!
                !Define step size
                t_step = minval([tau_bounce/20.d0,tau_collision/20.d0]);
!
                !Define number of steps
                n_time_steps = int(flight_time_multiplier &
                                   & * maxval([tau_collision,tau_bounce**2/tau_collision])/t_step)

                call write_transport_metadata_row(file_id_transport_metadata,nu_star,mag_axis_R0,aiota,q_saf, &
                                                   & vmod,coll_freq,t_step,n_time_steps,n_particles,seed_option, &
                                                   & n_field_periods)

                print *, 'Number of time steps (measurements)',n_time_steps
                print *, 'Total physical orbit flight time per particle',n_time_steps*t_step
!
                print *, 'Start calculation of diffusion coefficient'
!
                if(boole_write_particle_histories) then
                    call calc_flux_deviation(n_particles,n_time_steps,t_step,vmod,boole_collisions,boole_random_precalc, &
                                        & coll_freq,file_id_transp_diff_coef,file_id_psi2,file_id_std_psi2, &
                                        & i_integrator_type,idiffcoef_output,nu_star,file_id_loss_events, &
                                        & file_id_loss_summary,file_id_particle_histories)
                else
                    call calc_flux_deviation(n_particles,n_time_steps,t_step,vmod,boole_collisions,boole_random_precalc, &
                                        & coll_freq,file_id_transp_diff_coef,file_id_psi2,file_id_std_psi2, &
                                        & i_integrator_type,idiffcoef_output,nu_star,file_id_loss_events,file_id_loss_summary)
                endif
!
            enddo !n_nu_scans
            close(file_id_loss_events)
            close(file_id_loss_summary)
            close(file_id_transport_metadata)
            if(boole_write_particle_histories) close(file_id_particle_histories)
!
            select case(idiffcoef_output)
                case(1)
                    close(file_id_transp_diff_coef)
                case(2)
                    close(file_id_psi2)
                    close(file_id_std_psi2)
                case(3)
                    close(file_id_transp_diff_coef)
                    close(file_id_psi2)
                    close(file_id_std_psi2)
            end select
!
        end subroutine calc_mono_energetic_transp_coef_nu_scan
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine calc_numerical_diff_coef()
!
            use flux_deviation_mod, only: calc_flux_deviation
            use orbit_timestep_gorilla_mod, only: initialize_gorilla
            use constants, only: ev2erg
            use tetra_physics_mod, only: particle_mass
            use tetra_grid_settings_mod, only: n1,n2,n3,set_grid_size
            use fluxtv_mod, only: load_flux_tube_volume
            use mono_energetic_transp_coef_settings_mod, only: boole_collisions,boole_random_precalc,i_integrator_type, &
                                           & filename_delta_s_squared, filename_std_dvt_delta_s_squared, &
                                           & idiffcoef_output,total_MC_time,energy_eV, &
                                           & n_particles,filename_numerical_diff_coef, &
                                           & n_time_steps => nt_steps_numerical_diff
!
            implicit none
!
            integer                                     :: file_id_psi2,file_id_std_psi2,file_id_transp_diff_coef
            double precision                            :: coll_freq,t_step
!
            !Load fluxtube volume for a starting position (File to be chosen in gorilla_applets.inp)
            call set_grid_size([n1,n2,n3])      !Rectangular grid:                  [n1,n2,n3] = [nR,nphi,nZ]
                                                !Field-aligned grid:                [n1,n2,n3] = [ns,nphi,ntheta]
            call load_flux_tube_volume()
!
            !Initialize GORILLA
            call initialize_gorilla()
!
            !Compute velocity module from kinetic energy dependent on particle species
            vmod=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
!
            print *, 'Numerical radial transport diffusion coefficient:'
            print *, ''
            print *, 'Number of particles', n_particles
!
            select case(idiffcoef_output)
                case(1)
                    open(newunit=file_id_transp_diff_coef,file=filename_numerical_diff_coef)
                case(2)
                    open(newunit=file_id_psi2,file=filename_delta_s_squared)
                    open(newunit=file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
                case(3)
                    open(newunit=file_id_transp_diff_coef,file=filename_numerical_diff_coef)
                    open(newunit=file_id_psi2,file=filename_delta_s_squared)
                    open(newunit=file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
            end select
!
            !Define step size
            t_step = total_MC_time/dble(n_time_steps)
!
            !Collision frequency (not used)
            coll_freq = 0.d0
!
            print *, 'Number of time steps (measurements)',n_time_steps
!
            print *, 'Start calculation of diffusion coefficient'
!
            call calc_flux_deviation(n_particles,n_time_steps,t_step,vmod,boole_collisions,boole_random_precalc,coll_freq, &
                                & file_id_transp_diff_coef,file_id_psi2,file_id_std_psi2,i_integrator_type,idiffcoef_output)
!
            select case(idiffcoef_output)
                case(1)
                    close(file_id_transp_diff_coef)
                case(2)
                    close(file_id_psi2)
                    close(file_id_std_psi2)
                case(3)
                    close(file_id_transp_diff_coef)
                    close(file_id_psi2)
                    close(file_id_std_psi2)
            end select
!
        end subroutine calc_numerical_diff_coef
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        !Resolve q and iota on the marker surface. EFIT uses its symmetry-flux q
        !spline; VMEC uses the iota spline from the equilibrium. Must be called
        !after initialize_mono_energetic_transp_coef loads pos_fluxtv_mat.
        subroutine get_local_rotational_transform(q_saf,aiota)
!
            use tetra_physics_mod, only: coord_system
            use tetra_grid_settings_mod, only: grid_kind
            use fluxtv_mod, only: pos_fluxtv_mat
            use magdata_in_symfluxcoordinates_mod, only: magdata_in_symfluxcoord_ext
            use splint_vmec_data_mod, only: splint_vmec_data
            use boozer_chartmap_mod, only: chartmap_iota
!
            implicit none
!
            double precision, intent(out)               :: q_saf,aiota
            double precision                            :: s_start,theta_start,varphi_start
            !Outputs of magdata_in_symfluxcoord_ext that are not needed here
            double precision                            :: psi_pol,dq_ds,sqrtg,bmod,dbmod_dtheta
            double precision                            :: R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta
            !Outputs of splint_vmec_data that are not needed here
            double precision                            :: A_phi,A_theta,dA_phi_ds,dA_theta_ds,alam
            double precision                            :: dR_dt,dR_dp,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
!
            q_saf = 2.d0
            aiota = 0.5d0
!
            select case(grid_kind)
                case(2)
                    if(coord_system.ne.2) then
                        error stop 'grid_kind=2 transport requires coord_system=2 for local q'
                    endif
                    s_start = pos_fluxtv_mat(1,1)
                    theta_start = 0.d0
                    call magdata_in_symfluxcoord_ext(1,s_start,psi_pol,theta_start,q_saf,dq_ds, &
                                                 sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta, &
                                                 Z,dZ_ds,dZ_dtheta)
                    if(abs(q_saf) <= tiny(q_saf)) error stop 'local safety factor is zero'
                    aiota = 1.d0/q_saf
                case(3)
                    if(coord_system.ne.2) then
                        error stop 'grid_kind=3 transport requires coord_system=2 for local iota'
                    endif
                    s_start = pos_fluxtv_mat(1,1)
                    theta_start = pos_fluxtv_mat(1,2)
                    varphi_start = pos_fluxtv_mat(1,3)
                    call splint_vmec_data(s_start,theta_start,varphi_start,A_phi,A_theta, &
                                          & dA_phi_ds,dA_theta_ds,aiota,R,Z,alam,dR_ds, &
                                          & dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
                    if(abs(aiota) <= tiny(aiota)) error stop 'local rotational transform is zero'
                    q_saf = 1.d0/aiota
                case(6)
                    if(coord_system.ne.2) then
                        error stop 'grid_kind=6 transport requires coord_system=2 for local iota'
                    endif
                    s_start = pos_fluxtv_mat(1,1)
                    aiota = chartmap_iota(s_start)
                    if(abs(aiota) <= tiny(aiota)) error stop 'local rotational transform is zero'
                    q_saf = 1.d0/aiota
                case(1)
                    print *, 'WARNING: grid_kind = 1, cannot look up the local safety factor.'
                    print *, '         Falling back to the hardcoded q = 2 for the bounce and'
                    print *, '         collision times. Use grid_kind = 2 to use the real q(s).'
                case default
                    error stop 'transport normalization is undefined for this grid_kind'
            end select
!
        end subroutine get_local_rotational_transform
!
        subroutine write_transport_metadata_header(file_id,seed_filename)
!
            implicit none
!
            integer, intent(in) :: file_id
            character(*), intent(in) :: seed_filename
!
            write(file_id,'(a)') '# nu_star_standard = R0 * nu_collision / (abs(iota) * speed)'
            write(file_id,'(a,a)') '# random_seed_filename = ',trim(seed_filename)
            write(file_id,'(a)') '# nu_star R0_cm iota q speed_cm_s nu_collision_s-1 '// &
                                 & 'time_step_s n_steps n_particles seed_option n_field_periods'
!
        end subroutine write_transport_metadata_header
!
        subroutine write_transport_metadata_row(file_id,nu_star,major_radius,aiota,q_saf,speed, &
                                                & nu_collision,time_step,n_steps,n_particles,seed_option, &
                                                & n_field_periods)
!
            implicit none
!
            integer, intent(in) :: file_id,n_particles,seed_option,n_field_periods
            integer(kind=8), intent(in) :: n_steps
            double precision, intent(in) :: nu_star,major_radius,aiota,q_saf,speed
            double precision, intent(in) :: nu_collision,time_step
!
            write(file_id,*) nu_star,major_radius,aiota,q_saf,speed,nu_collision,time_step, &
                             & n_steps,n_particles,seed_option,n_field_periods
!
        end subroutine write_transport_metadata_row
!
    end module gorilla_applets_sub_mod
