# FortSym-generated pointwise kernels for the linear delta-f NTV estimator.
# Generator: gen/deltaf_kernels (not run at build time; output is committed).
add_library(GORILLA_APPLETS_DELTAF_KERNELS STATIC
  deltaf_dke_characteristics_kernel.f90
  deltaf_dke_birth_kernel.f90
  deltaf_dke_reference_kernel.f90
  deltaf_dke_local_kernel.f90
  deltaf_lerp_kernel.f90
  deltaf_dke_source_kernel.f90
  deltaf_dke_step_kernel.f90
  deltaf_fax_horner_kernel.f90
  deltaf_fax_energy_kernel.f90
  deltaf_fax_speed_kernel.f90
  deltaf_harmonic_kernel.f90
  deltaf_circular_chart_kernel.f90
  deltaf_weight_kernel.f90
  deltaf_fax_coordinates_kernel.f90
  deltaf_fax_chartmap_kernel.f90
  deltaf_fax_weight_kernel.f90
  )
set_target_properties(GORILLA_APPLETS_DELTAF_KERNELS PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(GORILLA_APPLETS_DELTAF_KERNELS PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR})
