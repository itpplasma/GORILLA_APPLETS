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

# No -ffpe-trap=invalid: LAPACK dgesv produces transient NaN in the
# make_tetra_physics OMP loop for grid_kind=5 with coord_system=2 (macOS ARM64).
add_compile_options(-g -fbacktrace -ffpe-trap=zero,overflow -fbounds-check -fopenmp -ffree-line-length-none)
add_link_options(-g -fbacktrace -ffpe-trap=zero,overflow -fbounds-check -fopenmp)

enable_testing()

add_subdirectory(SRC)
add_subdirectory(SRC_CORE)
add_subdirectory(TESTS)
add_subdirectory(SRC/generated)
add_subdirectory(TESTS/deltaf_kernels)

add_executable(gorilla_applets_main.x
  SRC/gorilla_applets_main.f90
)

target_link_libraries(gorilla_applets_main.x
  GORILLA_APPLETS
  )
target_include_directories(gorilla_applets_main.x PUBLIC ${CMAKE_BINARY_DIR}/SRC_CORE)
include(CTest)
if(BUILD_TESTING)
    add_executable(test_ou_collision_invariant.x
        SRC_CORE/constants_mod.f90
        SRC/collis_ions_mod.f90
        TESTS/test_ou_collision_invariant.f90)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ou_test_modules)
    set_target_properties(test_ou_collision_invariant.x PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ou_test_modules)
    target_include_directories(test_ou_collision_invariant.x PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/ou_test_modules)
    add_test(NAME ou_collision_invariant COMMAND test_ou_collision_invariant.x)
    # Analytic moment oracle; no equilibrium or orbit calculation is required.
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

if(BUILD_TESTING)
    add_executable(test_rmp_maxwellian_loading.x
        SRC/rmp_response_currents/rmp_maxwellian_loading_mod.f90
        TESTS/test_rmp_maxwellian_loading.f90)
    set_target_properties(test_rmp_maxwellian_loading.x PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/maxwellian_test_modules)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/maxwellian_test_modules)
    target_include_directories(test_rmp_maxwellian_loading.x PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/maxwellian_test_modules)
    add_test(NAME rmp_maxwellian_loading COMMAND test_rmp_maxwellian_loading.x)
endif()

if(BUILD_TESTING)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/potential_fixture)
    add_executable(test_rmp_equilibrium_potential.x
        TESTS/test_rmp_equilibrium_potential.f90)
    target_link_libraries(test_rmp_equilibrium_potential.x GORILLA_APPLETS)
    add_test(NAME rmp_equilibrium_potential COMMAND test_rmp_equilibrium_potential.x)
endif()

if(BUILD_TESTING)
    add_executable(test_canonical_ntv_kernel.x TESTS/test_canonical_ntv_kernel.f90)
    target_link_libraries(test_canonical_ntv_kernel.x GORILLA_APPLETS)
    add_test(NAME canonical_ntv_kernel COMMAND test_canonical_ntv_kernel.x)
    add_executable(test_ntv_lorentz_proposal.x TESTS/test_ntv_lorentz_proposal.f90)
    target_link_libraries(test_ntv_lorentz_proposal.x GORILLA_APPLETS)
    add_test(NAME ntv_lorentz_proposal COMMAND test_ntv_lorentz_proposal.x)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME canonical_ntv
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/TESTS/test_canonical_ntv.py)
    set_tests_properties(canonical_ntv PROPERTIES
        ENVIRONMENT
            "APPLETS_ROOT=${CMAKE_CURRENT_SOURCE_DIR};\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/canonical_ntv_work;\
OMP_NUM_THREADS=4"
        PROCESSORS 4
        TIMEOUT 600)
    add_executable(test_canonical_ntv_fax.x TESTS/test_canonical_ntv_fax.f90)
    target_link_libraries(test_canonical_ntv_fax.x GORILLA_APPLETS)
    add_test(NAME canonical_ntv_fax_kernel COMMAND test_canonical_ntv_fax.x)
    add_test(NAME canonical_ntv_fax
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/TESTS/test_canonical_ntv_fax.py)
    set_tests_properties(canonical_ntv_fax PROPERTIES
        ENVIRONMENT
            "APPLETS_ROOT=${CMAKE_CURRENT_SOURCE_DIR};\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/canonical_ntv_fax_work;\
OMP_NUM_THREADS=1"
        TIMEOUT 600)
    add_test(NAME canonical_ntv_lorentz
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/TESTS/test_canonical_ntv_lorentz.py)
    set_tests_properties(canonical_ntv_lorentz PROPERTIES
        ENVIRONMENT
            "APPLETS_ROOT=${CMAKE_CURRENT_SOURCE_DIR};\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/canonical_ntv_lorentz_work;\
OMP_NUM_THREADS=4"
        PROCESSORS 4
        TIMEOUT 600)
endif()

if(BUILD_TESTING)
    # Boozer chartmap of the analytic circular field (two-route NTV tests).
    add_executable(write_circular_boozer_chartmap.x
        TESTS/write_circular_boozer_chartmap.f90)
    target_link_libraries(write_circular_boozer_chartmap.x GORILLA netcdff
        netCDF::netcdf)
    target_include_directories(write_circular_boozer_chartmap.x PUBLIC
        ${CMAKE_BINARY_DIR}/SRC_CORE)
    add_executable(test_canonical_ntv_chartmap.x TESTS/test_canonical_ntv_chartmap.f90)
    target_link_libraries(test_canonical_ntv_chartmap.x GORILLA_APPLETS)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/canonical_ntv_chartmap_unit)
    add_test(NAME canonical_ntv_chartmap_unit
        COMMAND test_canonical_ntv_chartmap.x
                $<TARGET_FILE:write_circular_boozer_chartmap.x>
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/canonical_ntv_chartmap_unit)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME canonical_ntv_chartmap
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/TESTS/test_canonical_ntv_chartmap.py)
    set_tests_properties(canonical_ntv_chartmap PROPERTIES
        ENVIRONMENT
            "APPLETS_ROOT=${CMAKE_CURRENT_SOURCE_DIR};\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
CHARTMAP_WRITER_BIN=$<TARGET_FILE:write_circular_boozer_chartmap.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/canonical_ntv_chartmap_work;\
OMP_NUM_THREADS=4"
        PROCESSORS 4
        TIMEOUT 1200)
endif()

if(BUILD_TESTING)
    add_executable(test_rmp_gc_measure.x
        SRC/rmp_response_currents/rmp_gc_measure_mod.f90
        TESTS/test_rmp_gc_measure.f90)
    set_target_properties(test_rmp_gc_measure.x PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/measure_test_modules)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/measure_test_modules)
    target_include_directories(test_rmp_gc_measure.x PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/measure_test_modules)
    add_test(NAME rmp_gc_measure COMMAND test_rmp_gc_measure.x)
endif()

if(BUILD_TESTING)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME rmp_fm_stationarity
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/TESTS/test_rmp_fm_stationarity.py)
    set_tests_properties(rmp_fm_stationarity PROPERTIES
        ENVIRONMENT
            "APPLETS_ROOT=${CMAKE_CURRENT_SOURCE_DIR};\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/rmp_fm_stationarity_work;\
OMP_NUM_THREADS=4"
        PROCESSORS 4
        TIMEOUT 900)
endif()
