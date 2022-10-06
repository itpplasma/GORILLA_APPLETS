if ("x$ENV{GORILLA_PATH}" STREQUAL "x")
  set(GORILLA_PATH "${CMAKE_SOURCE_DIR}/../GORILLA/BUILD/")
else()
  set(GORILLA_PATH "$ENV{GORILLA_PATH}/")
endif()

find_library(GORILLA GORILLA ${GORILLA_PATH}/SRC/ NO_DEFAULT_PATH)

add_library(GORILLA_APPLETS STATIC
  alpha_lifetime_gorilla_mod.f90
  boltzmann_mod.f90
  direct_vmec_integrator_mod.f90
  flux_deviation_mod.f90
  fluxtv_mod.f90
  fluxtv_pusher_mod.f90
  gorilla_applets_settings_mod.f90
  gorilla_applets_sub_mod.f90
  llsq.f90
  mono_energetic_transp_coef_settings_mod.f90
  poincare_invariances_mod.f90
  sub_alpha_lifetime_can.f90
  total_dwell_times_mod.f90
  )
target_include_directories(GORILLA_APPLETS PUBLIC "${GORILLA_PATH}/SRC/")

add_executable(gorilla_applets_main.x
  gorilla_applets_main.f90
  )

target_link_libraries(gorilla_applets_main.x
  GORILLA_APPLETS
  ${GORILLA}
  netcdff
  netcdf
  lapack)
target_include_directories(gorilla_applets_main.x PUBLIC "${GORILLA_PATH}/SRC/")

set_target_properties (GORILLA_APPLETS PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(GORILLA_APPLETS PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
