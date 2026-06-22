cmake_minimum_required(VERSION 3.12)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (GORILLA_APPLETS
  VERSION 1.0.0)

enable_language(Fortran)

if(UNIX AND NOT APPLE)
  include_directories(/usr/include)
elseif(APPLE)
  # Locate netCDF via nf-config / nc-config so the build works with either
  # Homebrew or MacPorts instead of assuming a fixed prefix. netCDF-Fortran
  # (libnetcdff, the .mod) and netCDF-C (libnetcdf) live under different
  # prefixes on Homebrew, so query both. Fall back to the common prefixes
  # when the config tools are not on PATH.
  find_program(NF_CONFIG nf-config)
  find_program(NC_CONFIG nc-config)
  if(NF_CONFIG)
    execute_process(COMMAND ${NF_CONFIG} --includedir
                    OUTPUT_VARIABLE NETCDF_FORTRAN_INCLUDEDIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${NF_CONFIG} --prefix
                    OUTPUT_VARIABLE NETCDF_FORTRAN_PREFIX
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    include_directories(${NETCDF_FORTRAN_INCLUDEDIR})
    link_directories(${NETCDF_FORTRAN_PREFIX}/lib)
  endif()
  if(NC_CONFIG)
    execute_process(COMMAND ${NC_CONFIG} --includedir
                    OUTPUT_VARIABLE NETCDF_C_INCLUDEDIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${NC_CONFIG} --libdir
                    OUTPUT_VARIABLE NETCDF_C_LIBDIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    include_directories(${NETCDF_C_INCLUDEDIR})
    link_directories(${NETCDF_C_LIBDIR})
  endif()
  if(NOT NF_CONFIG AND NOT NC_CONFIG)
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
      include_directories(/opt/homebrew/include /opt/local/include)
      link_directories(/opt/homebrew/lib /opt/local/lib)
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
      include_directories(/usr/local/include)
      link_directories(/usr/local/lib)
    endif()
  endif()
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(netCDF REQUIRED)

# netCDF-Fortran ships a CMake config in upstream installs (AFS modules,
# Homebrew, MacPorts) but not in Debian's libnetcdff-dev. Fall back to a plain
# library lookup so SRC/CMakeLists.txt can still link `netCDF::netcdff`.
find_package(netCDF-Fortran QUIET CONFIG)
if(NOT netCDF-Fortran_FOUND)
  find_library(NETCDFF_LIBRARY NAMES netcdff REQUIRED)
  add_library(netCDF::netcdff UNKNOWN IMPORTED)
  set_target_properties(netCDF::netcdff PROPERTIES
    IMPORTED_LOCATION "${NETCDFF_LIBRARY}")
endif()

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