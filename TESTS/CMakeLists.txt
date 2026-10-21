add_executable(test_multi_species_collision.x
    test_multi_species_collision.f90
)
target_link_libraries(test_multi_species_collision.x GORILLA_APPLETS)
target_include_directories(test_multi_species_collision.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)

add_test(NAME multi_species_collision COMMAND test_multi_species_collision.x)

add_executable(test_transport_loss_moments.x
    test_transport_loss_moments.f90
)
target_link_libraries(test_transport_loss_moments.x GORILLA_APPLETS)
target_include_directories(test_transport_loss_moments.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)

add_test(NAME transport_loss_moments COMMAND test_transport_loss_moments.x)

add_executable(test_transport_normalization.x
    test_transport_normalization.f90
)
target_link_libraries(test_transport_normalization.x GORILLA_APPLETS)
target_include_directories(test_transport_normalization.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)

add_test(NAME transport_normalization COMMAND test_transport_normalization.x)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(
    NAME alpha_lifetime
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_alpha_lifetime.py
)
set_tests_properties(alpha_lifetime PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/alpha_lifetime_work"
)

add_test(
    NAME field_line_tracing
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_field_line_tracing.py
)
set_tests_properties(field_line_tracing PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/field_line_tracing_work"
)

add_test(
    NAME mono_energetic_transport_slow
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_mono_energetic_transport_slow.py
)
set_tests_properties(mono_energetic_transport_slow PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/mono_energetic_transport_work"
)

add_test(
    NAME helical_core
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_helical_core.py
)
set_tests_properties(helical_core PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/helical_core_work"
)

add_test(
    NAME helical_core_rectangular
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_helical_core.py
)
set_tests_properties(helical_core_rectangular PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/helical_core_rectangular_work;\
HELICAL_CORE_GRID_KIND=1"
)

add_executable(write_chartmap_fixture.x
    write_chartmap_fixture.f90
)
target_link_libraries(write_chartmap_fixture.x netcdf netcdff)

add_test(
    NAME mono_energetic_transport_chartmap
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_mono_energetic_transport_chartmap.py
)
set_tests_properties(mono_energetic_transport_chartmap PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
CHARTMAP_FIXTURE_BIN=$<TARGET_FILE:write_chartmap_fixture.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/mono_energetic_transport_chartmap_work"
)

add_executable(test_flux_mesh_potential_driver.x
    test_flux_mesh_potential_driver.f90
)
target_link_libraries(test_flux_mesh_potential_driver.x GORILLA_APPLETS)
target_include_directories(test_flux_mesh_potential_driver.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)

add_test(
    NAME flux_mesh_potential_chartmap
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_flux_mesh_potential.py
)
set_tests_properties(flux_mesh_potential_chartmap PROPERTIES
    ENVIRONMENT
        "GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
DRIVER_BIN=$<TARGET_FILE:test_flux_mesh_potential_driver.x>;\
CHARTMAP_WRITER_BIN=$<TARGET_FILE:write_circular_boozer_chartmap.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/flux_mesh_potential_work"
)

add_executable(test_local_collision_background.x test_local_collision_background.f90)
target_link_libraries(test_local_collision_background.x GORILLA_APPLETS)
add_test(NAME local_collision_background
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_local_collision_background.py
        $<TARGET_FILE:test_local_collision_background.x>
        $<TARGET_FILE:write_circular_boozer_chartmap.x>
        ${CMAKE_CURRENT_BINARY_DIR}/local_collision_background_work)
set_tests_properties(local_collision_background PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1"
    PROCESSORS 2 TIMEOUT 120)

add_test(NAME deltaf_3d_orbits
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_deltaf_3d.py)
set_tests_properties(deltaf_3d_orbits PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;CHARTMAP_WRITER_BIN=$<TARGET_FILE:write_circular_boozer_chartmap.x>;WORK_DIR=${CMAKE_BINARY_DIR}/deltaf_3d_work;OMP_NUM_THREADS=1"
    PROCESSORS 4 TIMEOUT 600)

add_executable(deltaf_3d_orbit_probe.x deltaf_3d_orbit_probe.f90)
target_link_libraries(deltaf_3d_orbit_probe.x GORILLA_APPLETS)

add_test(NAME deltaf_3d_ripple_wells
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_deltaf_3d_wells.py)
set_tests_properties(deltaf_3d_ripple_wells PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;CHARTMAP_WRITER_BIN=$<TARGET_FILE:write_circular_boozer_chartmap.x>;ORBIT_PROBE_BIN=$<TARGET_FILE:deltaf_3d_orbit_probe.x>;WORK_DIR=${CMAKE_BINARY_DIR}/deltaf_3d_wells;OMP_NUM_THREADS=1"
    TIMEOUT 300)

add_test(NAME deltaf_3d_linear_limit
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_deltaf_3d_linear.py)
set_tests_properties(deltaf_3d_linear_limit PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;CHARTMAP_WRITER_BIN=$<TARGET_FILE:write_circular_boozer_chartmap.x>;WORK_DIR=${CMAKE_BINARY_DIR}/deltaf_3d_linear;OMP_NUM_THREADS=1"
    TIMEOUT 300)
add_test(NAME canonical_species_density
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_canonical_species_density.py)
set_tests_properties(canonical_species_density PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/canonical_species_density_work;\
OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1"
    TIMEOUT 120)

# Link the reader probe to the same core library as the production executable.
add_executable(test_marker_species.x ../SRC_CORE/TESTS/test_marker_species.f90)
target_link_libraries(test_marker_species.x GORILLA_APPLETS)
target_include_directories(test_marker_species.x PRIVATE ${CMAKE_BINARY_DIR}/SRC_CORE)
add_test(NAME marker_species
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/../SRC_CORE/TESTS/test_marker_species.py
        $<TARGET_FILE:test_marker_species.x>
        ${CMAKE_CURRENT_SOURCE_DIR}/../SRC_CORE/../INPUT/gorilla.inp)
