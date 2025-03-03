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
                                    & calc_numerical_diff_coef,calc_mono_energetic_transp_coef_nu_scan
!
    contains
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
               ds_dr = 2.d0*sqrt(s)*sqrt(pi*bmod00/torflux)
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
                                       & idiffcoef_output, energy_eV, v_E, nu_star, n_particles
!
            implicit none
!
            integer                                     :: file_id_psi2,file_id_std_psi2,file_id_transp_diff_coef
            integer(kind=8)                             :: n_time_steps
            double precision                            :: tau_bounce,coll_freq,tau_collision,t_step
!
            !Initialize GORILLA with electrostatic potential in accordance with Mach number for mono-energetic transport coefficient
            call initialize_mono_energetic_transp_coef()

            !File ids for squared deviaton, standard deviation and diffusion coefficient
            file_id_psi2 = 14
            file_id_std_psi2 = 15
            file_id_transp_diff_coef = 16
!
            !nu_star = [(2.d0**i, i=nu_start,(nu_start-(n_nu_scans+1)), -1)] ! nu_star = \frac{R_0 \nu_c}{\iota v_{mod}}
!
            select case(idiffcoef_output)
                case(1)
                    open(file_id_transp_diff_coef,file=filename_transp_diff_coef)
                case(2)
                    open(file_id_psi2,file=filename_delta_s_squared)
                    open(file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
                case(3)
                    open(file_id_transp_diff_coef,file=filename_transp_diff_coef)
                    open(file_id_psi2,file=filename_delta_s_squared)
                    open(file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
            end select
!
            print *, 'Mono-energetic radial transport coefficient:'
            print *, ''
            print *, 'Normalized collisionality', nu_star
            print *, 'Mach number', v_E
            print *, 'Number of particles', n_particles
!
            select case(grid_kind)
                case(1,2)
                    !Define bounce time for Tokamak
                    tau_bounce = 2.d0*pi*mag_axis_R0/vmod*2.d0  !q=2
                    tau_collision = mag_axis_R0*2.d0/(vmod*nu_star)
                case(3)
                    !Define bounce time for Stellarator
                    tau_bounce = 2.d0*pi*mag_axis_R0/vmod/n_field_periods
                    tau_collision = mag_axis_R0/(vmod*nu_star)
            end select
!
            !Define collision frequency from collision time
            coll_freq = 1.d0/tau_collision
!
            !Define step size
            t_step = minval([tau_bounce/20.d0,tau_collision/20.d0]);
!
            !Define number of steps
            n_time_steps = int(10.d0* maxval([tau_collision,tau_bounce**2/tau_collision])/t_step) !(not couting initial measurement)

            print *, 'Number of time steps (measurements)',n_time_steps
!
            print *, 'Start calculation of diffusion coefficient'
!
            call calc_flux_deviation(n_particles,n_time_steps,t_step,vmod,boole_collisions,boole_random_precalc,coll_freq, &
                                & file_id_transp_diff_coef,file_id_psi2,file_id_std_psi2,i_integrator_type,idiffcoef_output, &
                                & nu_star)
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
                                       & nu_exp_basis
!
            implicit none
!
            integer                                     :: i,file_id_psi2,file_id_std_psi2,file_id_transp_diff_coef
            integer(kind=8)                             :: n_time_steps
            double precision                            :: tau_bounce,coll_freq,tau_collision,t_step, nu_star
!
            !Initialize GORILLA with electrostatic potential in accordance with Mach number for mono-energetic transport coefficient
            call initialize_mono_energetic_transp_coef()

            !File ids for squared deviaton, standard deviation and diffusion coefficient
            file_id_psi2 = 14
            file_id_std_psi2 = 15
            file_id_transp_diff_coef = 16
!
            select case(idiffcoef_output)
                case(1)
                    open(file_id_transp_diff_coef,file=filename_transp_diff_coef)
                case(2)
                    open(file_id_psi2,file=filename_delta_s_squared)
                    open(file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
                case(3)
                    open(file_id_transp_diff_coef,file=filename_transp_diff_coef)
                    open(file_id_psi2,file=filename_delta_s_squared)
                    open(file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
            end select
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
                    case(1,2)
                        !Define bounce time for Tokamak
                        tau_bounce = 2.d0*pi*mag_axis_R0/vmod*2.d0  !q=2
                        tau_collision = mag_axis_R0*2.d0/(vmod*nu_star)
                    case(3)
                        !Define bounce time for Stellarator
                        tau_bounce = 2.d0*pi*mag_axis_R0/vmod/n_field_periods
                        tau_collision = mag_axis_R0/(vmod*nu_star)
                end select
!
                !Define collision frequency from collision time
                coll_freq = 1.d0/tau_collision
!
                !Define step size
                t_step = minval([tau_bounce/20.d0,tau_collision/20.d0]);
!
                !Define number of steps
                n_time_steps = int(10.d0* maxval([tau_collision,tau_bounce**2/tau_collision])/t_step) !(not couting initial measurement)

                print *, 'Number of time steps (measurements)',n_time_steps
                print *, 'Total physical orbit flight time per particle',n_time_steps*t_step
!
                print *, 'Start calculation of diffusion coefficient'
!
                call calc_flux_deviation(n_particles,n_time_steps,t_step,vmod,boole_collisions,boole_random_precalc,coll_freq, &
                                    & file_id_transp_diff_coef,file_id_psi2,file_id_std_psi2,i_integrator_type,idiffcoef_output, &
                                    & nu_star)
!
            enddo !n_nu_scans
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
            !file id for squared deviaton
            file_id_psi2 = 14
            file_id_std_psi2 = 15
            file_id_transp_diff_coef = 16
!
            print *, 'Numerical radial transport diffusion coefficient:'
            print *, ''
            print *, 'Number of particles', n_particles
!
            select case(idiffcoef_output)
                case(1)
                    open(file_id_transp_diff_coef,file=filename_numerical_diff_coef)
                case(2)
                    open(file_id_psi2,file=filename_delta_s_squared)
                    open(file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
                case(3)
                    open(file_id_transp_diff_coef,file=filename_numerical_diff_coef)
                    open(file_id_psi2,file=filename_delta_s_squared)
                    open(file_id_std_psi2,file=filename_std_dvt_delta_s_squared)
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
    end module gorilla_applets_sub_mod
