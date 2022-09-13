!
    module gorilla_applets_settings_mod
!
        implicit none
!
        private
!
        !Fluxtube volume OR mono-energetic radial transport coefficient OR numerical diffusion coefficient
        integer, public, protected :: i_option
!
        !Options for pre-computation of fluxtube volume
        character(1024),public,protected :: filename_fluxtv_precomp
        double precision, public, protected :: start_pos_x1,start_pos_x2,start_pos_x3,t_step_fluxtv,energy_eV_fluxtv
        integer(kind=8), public, protected :: nt_steps_fluxtv
!
        !Settings of simulation for transport coefficients
        integer, public, protected :: i_integrator_type, n_particles
        logical, public, protected :: boole_collisions,boole_psi_mat
        double precision, public, protected :: energy_eV
        character(1024),public,protected :: filename_fluxtv_load
!
        !Settings for random numbers to model collisions
        logical, public, protected :: boole_random_precalc
        integer, public, protected :: seed_option
!
        !Settings for output of diffusion coefficient and standard deviation of diffusion coefficient
        integer, public, protected :: idiffcoef_output
        character(50),public,protected :: filename_transp_diff_coef, filename_delta_s_squared, filename_std_dvt_delta_s_squared
!
        !Settings for mono-energetic transport coefficient for given collisionality and Mach number
        character(50),public,protected :: filename_numerical_diff_coef
        double precision, public, protected :: v_E,nu_star
!
        !Settings for scan over collisionality: Mono-energetic transport coefficent for constant Mach number
        double precision, public, protected :: nu_star_start,nu_exp_basis
        integer, public, protected :: n_nu_scans
!
        !Settings for numerical diffusion coefficient (no collisions)
        double precision, public, protected :: lambda_start, total_MC_time
        integer(kind=8), public, protected :: nt_steps_numerical_diff
!
        public :: load_gorilla_applets_inp

        NAMELIST /GORILLA_APPLETS_NML/ i_option, filename_fluxtv_precomp, start_pos_x1, start_pos_x2, start_pos_x3, t_step_fluxtv, &
            & i_integrator_type, n_particles, boole_collisions, energy_eV, boole_random_precalc, seed_option, idiffcoef_output, &
            & filename_transp_diff_coef, filename_delta_s_squared, filename_std_dvt_delta_s_squared, v_E, nu_star, lambda_start, &
            & total_MC_time, energy_eV_fluxtv,filename_numerical_diff_coef,filename_fluxtv_load,nt_steps_numerical_diff, &
            & nt_steps_fluxtv,boole_psi_mat,nu_star_start,nu_exp_basis,n_nu_scans
!
    contains
!            
        subroutine load_gorilla_applets_inp()
!
            open(unit=11, file='gorilla_applets.inp', status='unknown')
            read(11,nml=GORILLA_APPLETS_NML)
            close(11)

            print *,'Mono-energetic transport coefficient: Loaded input data from input file'
!            
        end subroutine load_gorilla_applets_inp

    end module gorilla_applets_settings_mod
!    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
