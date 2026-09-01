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
        character(256), public, protected :: random_seed_filename = 'seed.inp'
!
        !Settings for output of diffusion coefficient and standard deviation of diffusion coefficient
        integer, public, protected :: idiffcoef_output
        character(50),public,protected :: filename_transp_diff_coef, filename_delta_s_squared, filename_std_dvt_delta_s_squared
!
        !Settings for mono-energetic transport coefficient for given collisionality and Mach number
        character(50),public,protected :: filename_numerical_diff_coef
        double precision, public, protected :: v_E,nu_star
        double precision, public, protected :: flight_time_multiplier = 10.d0
        logical, public, protected :: boole_write_particle_histories = .false.
        character(256), public, protected :: filename_particle_histories = 'particle_histories.dat'
        character(256), public, protected :: filename_transport_metadata = 'transport_metadata.dat'
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
            & boole_psi_mat,nu_star_start,nu_exp_basis,n_nu_scans,flight_time_multiplier, &
            & random_seed_filename,boole_write_particle_histories,filename_particle_histories, &
            & filename_transport_metadata
!
    contains
!            
        subroutine load_mono_energetic_transp_coef_inp()
!
            integer :: inp_unit
!
            open(newunit=inp_unit, file='mono_energetic_transp_coef.inp', status='unknown')
            read(inp_unit,nml=TRANSPCOEFNML)
            close(inp_unit)

            if(flight_time_multiplier <= 0.d0) then
                error stop 'flight_time_multiplier must be positive'
            endif

            if(boole_write_particle_histories .and. .not.boole_psi_mat) then
                error stop 'boole_write_particle_histories requires boole_psi_mat'
            endif

            print *,'Mono-energetic transport coefficient: Loaded input data from input file'
!            
        end subroutine load_mono_energetic_transp_coef_inp

    end module mono_energetic_transp_coef_settings_mod
!    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
