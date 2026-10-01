cmake_minimum_required(VERSION 3.12)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (GORILLA_APPLETS
  VERSION 1.0.0)

enable_language(Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(netCDF REQUIRED)

# netCDF-Fortran shim: on Debian/Ubuntu, libnetcdff-dev ships the library and
# .mod files but NOT netCDF-FortranConfig.cmake (upstream only generates it for
# CMake-based builds; Debian builds with autotools). Try find_package first so
# machines with a proper config file (source builds, conda, other distros) keep
# using it; otherwise synthesize the netCDF::netcdff IMPORTED target manually
# from the discovered .so and .mod paths.
find_package(netCDF-Fortran QUIET CONFIG)
if(NOT netCDF-Fortran_FOUND)
    message(STATUS "netCDF-Fortran: no CMake config file found, using manual discovery")
    find_library(NETCDFF_LIBRARY NAMES netcdff REQUIRED)
    find_path(NETCDFF_MODULE_DIR NAMES netcdf.mod
              HINTS /usr/include /usr/local/include
              PATH_SUFFIXES netcdf-fortran
              REQUIRED)
    add_library(netCDF::netcdff UNKNOWN IMPORTED)
    set_target_properties(netCDF::netcdff PROPERTIES
        IMPORTED_LOCATION "${NETCDFF_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${NETCDFF_MODULE_DIR}")
    message(STATUS "netCDF-Fortran library : ${NETCDFF_LIBRARY}")
    message(STATUS "netCDF-Fortran modules : ${NETCDFF_MODULE_DIR}")
endif()

get_target_property(_NETCDFF_INC netCDF::netcdff INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${netCDF_INCLUDE_DIR} ${_NETCDFF_INC})

add_compile_options(-g -fbacktrace -ffpe-trap=zero,overflow -fbounds-check -fopenmp)
add_link_options(-g -fbacktrace -ffpe-trap=zero,overflow -fbounds-check -fopenmp)
# NOTE: For coord_system=2 (SFL) with grid_kind=5 (analytic circular tokamak),
# gfortran's -ffpe-trap=invalid triggers on transient NaN inside LAPACK dgesv
# during the make_tetra_physics OMP loop. This is a known gfortran+OMP+LAPACK
# interaction on macOS ARM64 and does not indicate a physics bug (the analytic-
# circ SFL orbit check reproduces q_analytic to 1.14%, identical to coord_system=1).
# If using -ffpe-trap=invalid causes a SIGILL at !$OMP PARALLEL in
# calc_volume_integrals, rebuild with -ffpe-trap=zero,overflow only.

add_subdirectory(SRC)
add_subdirectory(SRC_CORE)

add_executable(gorilla_applets_main.x
  SRC/gorilla_applets_main.f90
)

target_link_libraries(gorilla_applets_main.x
  GORILLA_APPLETS
  )
target_include_directories(gorilla_applets_main.x PUBLIC ${CMAKE_BINARY_DIR}/SRC_CORE)
# Analytic moment oracle; no equilibrium or orbit calculation is required.
include(CTest)
if(BUILD_TESTING)
    add_executable(test_rmp_profile_moments.x
        SRC/rmp_response_currents/rmp_profile_moments_mod.f90
        TESTS/test_rmp_profile_moments.f90)
    set_target_properties(test_rmp_profile_moments.x PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/moment_test_modules)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/moment_test_modules)
    target_include_directories(test_rmp_profile_moments.x PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/moment_test_modules)
    add_test(NAME rmp_profile_moments COMMAND test_rmp_profile_moments.x)
endif()

if(BUILD_TESTING)
    add_executable(test_rmp_volume_loading.x
        SRC/rmp_response_currents/rmp_volume_loading_mod.f90
        TESTS/test_rmp_volume_loading.f90)
    set_target_properties(test_rmp_volume_loading.x PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/volume_test_modules)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/volume_test_modules)
    target_include_directories(test_rmp_volume_loading.x PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/volume_test_modules)
    add_test(NAME rmp_volume_loading COMMAND test_rmp_volume_loading.x)
endif()
