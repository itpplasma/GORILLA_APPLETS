!
  program gorilla_applets_main
!
  use gorilla_settings_mod, only: load_gorilla_inp,boole_helical_pert
  use tetra_grid_settings_mod, only: load_tetra_grid_inp
  use fluxtv_mod, only: calc_flux_tube_volume
  use gorilla_applets_sub_mod, only: calc_mono_energetic_transp_coef, &
                                    & calc_numerical_diff_coef, calc_mono_energetic_transp_coef_nu_scan
  use gorilla_applets_settings_mod, only: i_option,load_gorilla_applets_inp
  use alpha_lifetime_gorilla_mod, only: calc_alpha_lifetime_gorilla
  use direct_vmec_integrator_mod, only: direct_vmec_integrator
  use poincare_invariances_mod, only: compute_first_poincare_invariance
  use test_grid_mod, only: test_grid,test_west,test_eirene_grid
  use total_dwell_times_mod, only: calc_total_dwell_times
  use boltzmann_mod, only: calc_boltzmann
!
  implicit none
!
!-------------------------------------------------------------------------------------------!
!------------------------------------ Load input files  ------------------------------------!
!
    print *, ''
    print *, 'GORILLA - Guiding-center ORbit Integration with Local Linearization Approach'
    print *, 'Monte Carlo evaluation of mono-energetic transport coefficients'
    print *, ''
!
    !Load tetrahedronal grid input data
    call load_tetra_grid_inp()
!
    !Load GORILLA settings input data
    call load_gorilla_inp()
!
    !Load input from gorilla_applets.inp (v_E, n_particles)
    call load_gorilla_applets_inp()
!
!-------------------------------------------------------------------------------------------!
!-------------- Monte Carlo evaluation of mono-energetic transport coefficient -------------!
!
    !Options:
    !1 ... Pre-computation of fluxtube volume for GORILLA
    !2 ... Computation of mono-energetic transport coefficient - ONE given collisionality
    !3 ... Computation of mono-energetic transport coefficient - Scan over collisionalities
    !4 ... Computation of numerical diffusion coefficient (No collisions)
!
    select case(i_option)
!
!-------------------------------------------------------------------------------------------!
!
        case(1) !Pre-computation of fluxtube volume
                ! - GORILLA uses special piecewise linear elecromagnetic fields
                ! - Precompuation of fluxtube volume with GORILLA is needed for Monte Carlo
                !   evaluation of transport coefficients with GORILLA
!
            !No analytical perturbation of Tokamak equilibrium
            boole_helical_pert = .false.
!
            call calc_flux_tube_volume()
!
!-------------------------------------------------------------------------------------------!
!
        case(2) !Computation of mono-energetic transport coefficient - ONE given collisionality
!
            call calc_mono_energetic_transp_coef()
!
!-------------------------------------------------------------------------------------------!
!
        case(3) !Computation of mono-energetic transport coefficient - Scan over collisionalities
!
            call calc_mono_energetic_transp_coef_nu_scan()
!
!-------------------------------------------------------------------------------------------!
!
        case(4) !Computation of numerical diffusion coefficient (No collisions)
!
            call calc_numerical_diff_coef()
!
!-------------------------------------------------------------------------------------------!
!
        case(5) !Alpha Lifetime
!
            call calc_alpha_lifetime_gorilla
!
!-------------------------------------------------------------------------------------------!
!
        case(6) !DIRECT VMEC Integrator
!
            call direct_vmec_integrator
!
!-------------------------------------------------------------------------------------------!
!
        case(7) !Compute 1st Poincaré invariance
!
            call compute_first_poincare_invariance
!
!-------------------------------------------------------------------------------------------!
!
        case(8) !Test west
!
            call test_eirene_grid()
!
!-------------------------------------------------------------------------------------------!
!
        case(9) !Compute total dwell times
                !
            call calc_total_dwell_times
!
!-------------------------------------------------------------------------------------------!
!
        case(10) !Do the boltzmanntest
                !
            call calc_boltzmann
!
!-------------------------------------------------------------------------------------------!
!
    end select
!
  end program gorilla_applets_main
