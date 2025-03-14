!
  program gorilla_applets_main
!
  use gorilla_settings_mod, only: load_gorilla_inp,boole_helical_pert
  use tetra_grid_settings_mod, only: load_tetra_grid_inp
  use fluxtv_mod, only: calc_flux_tube_volume
  use gorilla_applets_sub_mod, only: calc_mono_energetic_transp_coef, &
                                    & calc_numerical_diff_coef, calc_mono_energetic_transp_coef_nu_scan
  use gorilla_applets_settings_mod, only: i_option,load_gorilla_applets_inp
  use mono_energetic_transp_coef_settings_mod, only: load_mono_energetic_transp_coef_inp
  use alpha_lifetime_gorilla_mod, only: calc_alpha_lifetime_gorilla
  use direct_vmec_integrator_mod, only: direct_vmec_integrator
  use poincare_invariances_mod, only: compute_first_poincare_invariance
  use reversibility_test_mod, only: make_reversibility_test
  use boltzmann_mod, only: calc_boltzmann
  use field_line_tracing_mod, only: calc_field_lines
  use divertor_heat_loads_mod, only: calc_divertor_heat_loads
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
    !Load input from gorilla_applets.inp (i_option)
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
!
            !Load input from mono_energetic_transp_coef.inp (v_E, n_particles)
            call load_mono_energetic_transp_coef_inp()
!
            call calc_mono_energetic_transp_coef()
!
!-------------------------------------------------------------------------------------------!
!
        case(3) !Computation of mono-energetic transport coefficient - Scan over collisionalities
!
!
            !Load input from mono_energetic_transp_coef.inp (v_E, n_particles)
            call load_mono_energetic_transp_coef_inp()
!
            call calc_mono_energetic_transp_coef_nu_scan()
!
!-------------------------------------------------------------------------------------------!
!
        case(4) !Computation of numerical diffusion coefficient (No collisions)
!
!
            !Load input from mono_energetic_transp_coef.inp (v_E, n_particles)
            call load_mono_energetic_transp_coef_inp()
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
        case(8) !Perform a reversibility test (demonstrate negative time)
!
            call make_reversibility_test()
!
!-------------------------------------------------------------------------------------------!
!
        case(9) !Do the boltzmanntest
!
            call calc_boltzmann
!
!-------------------------------------------------------------------------------------------!
!
        case(10) !Do field line tracing
!
            call calc_field_lines
!
!-------------------------------------------------------------------------------------------!
!
        case(11) !Calculate divertor heat loads
!
            call calc_divertor_heat_loads
!
!-------------------------------------------------------------------------------------------!
!
            
    end select
!
  end program gorilla_applets_main
