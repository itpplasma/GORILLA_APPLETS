add_library(GORILLA_APPLETS STATIC
  alpha_lifetime_gorilla_mod.f90
  boltzmann_mod.f90
  boltzmann_types_mod.f90
  collis_ions_mod.f90
  direct_vmec_integrator_mod.f90
  divertor_heat_loads_mod.f90
  field_line_tracing_mod.f90
  flux_deviation_mod.f90
  fluxtv_mod.f90
  fluxtv_pusher_mod.f90
  gorilla_applets_settings_mod.f90
  gorilla_applets_sub_mod.f90
  llsq.f90
  mono_energetic_transp_coef_settings_mod.f90
  poincare_invariances_mod.f90
  pusher_tetra_field_lines_mod.f90
  reversibility_test_mod.f90
  sub_alpha_lifetime_can.f90
  utils_data_pre_and_post_processing_mod.f90
  utils_orbit_timestep_mod.f90
  utils_parallelised_particle_pushing_mod.f90
  utils_write_data_to_files_mod.f90
  volume_integrals_and_sqrt_g_mod.f90
  zzg.f90
  )
target_link_libraries(GORILLA_APPLETS GORILLA)
target_include_directories(GORILLA_APPLETS PUBLIC ${CMAKE_BINARY_DIR}/SRC_CORE)

set_target_properties (GORILLA_APPLETS PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(GORILLA_APPLETS PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
