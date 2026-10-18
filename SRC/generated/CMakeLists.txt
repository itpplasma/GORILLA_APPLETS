# FortSym-generated pointwise kernels for the linear delta-f NTV estimator.
# Generator: gen/deltaf_kernels (not run at build time; output is committed).
add_library(GORILLA_APPLETS_DELTAF_KERNELS STATIC
  deltaf_harmonic_kernel.f90
  deltaf_circular_chart_kernel.f90
  deltaf_circular_flux_label_kernel.f90
  deltaf_canonical_momentum_kernel.f90
  deltaf_weight_kernel.f90
  )
set_target_properties(GORILLA_APPLETS_DELTAF_KERNELS PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(GORILLA_APPLETS_DELTAF_KERNELS PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR})
