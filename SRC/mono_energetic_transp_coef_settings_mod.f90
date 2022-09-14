!
    module mono_energetic_transp_coef_settings_mod
!
        implicit none
!
        private
!
        !Settings of simulation for transport coefficients
        integer, public, protected :: i_integrator_type, n_particles
        logical, public, protected :: boole_collisions,boole_psi_mat
        double precision, public, protected :: energy_eV
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
        public :: load_mono_energetic_transp_coef_inp

        NAMELIST /TRANSPCOEFNML/ &
            & i_integrator_type, n_particles, boole_collisions, energy_eV, boole_random_precalc, seed_option, idiffcoef_output, &
            & filename_transp_diff_coef, filename_delta_s_squared, filename_std_dvt_delta_s_squared, v_E, nu_star, lambda_start, &
            & total_MC_time,filename_numerical_diff_coef,nt_steps_numerical_diff, &
            & boole_psi_mat,nu_star_start,nu_exp_basis,n_nu_scans
!
    contains
!            
        subroutine load_mono_energetic_transp_coef_inp()
!
            open(unit=11, file='mono_energetic_transp_coef.inp', status='unknown')
            read(11,nml=TRANSPCOEFNML)
            close(11)

            print *,'Mono-energetic transport coefficient: Loaded input data from input file'
!            
        end subroutine load_mono_energetic_transp_coef_inp

    end module mono_energetic_transp_coef_settings_mod
!    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
