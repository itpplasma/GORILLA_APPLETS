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

foreach(name IN ITEMS homma_coefficients homma_kick stost_cartesian)
    add_executable(test_${name}.x test_${name}.f90)
    target_link_libraries(test_${name}.x GORILLA_APPLETS)
    target_include_directories(test_${name}.x
        PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)
    add_test(NAME ${name}
             COMMAND test_${name}.x ${CMAKE_CURRENT_SOURCE_DIR}/REFERENCE/homma)
endforeach()

find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Only where Wolfram Language is installed (not in CI); about 6 s. The proofs
# (DERIVATIONS/homma_hatayama_proofs.wl, about 3 min) are run by hand.
find_program(WOLFRAM_MATH math)
if(WOLFRAM_MATH)
    add_test(
        NAME homma_derivation_current
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/test_homma_derivation_current.py
    )
    set_tests_properties(homma_derivation_current PROPERTIES
        TIMEOUT 60
        ENVIRONMENT "APPLETS_ROOT=${CMAKE_SOURCE_DIR};WOLFRAM_MATH=${WOLFRAM_MATH}"
    )
endif()

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

add_test(
    NAME helical_core_conserving_collisions
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_helical_core_conserving.py
)
set_tests_properties(helical_core_conserving_collisions PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/helical_core_conserving_work"
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

add_test(
    NAME self_consistent_electric_field
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_self_consistent_electric_field.py
)
set_tests_properties(self_consistent_electric_field PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/self_consistent_electric_field_work"
)

# Reference tests of the applets without other end-to-end coverage
foreach(mode IN ITEMS boltzmann divertor_heat_loads anomalous_transport)
    add_test(
        NAME reference_${mode}
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/test_applet_reference.py ${mode}
    )
    set_tests_properties(reference_${mode} PROPERTIES
        ENVIRONMENT
            "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_ROOT=${CMAKE_SOURCE_DIR}/../GORILLA;\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/reference_${mode}_work"
    )
endforeach()

add_test(
    NAME thermal_force_velocity
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_thermal_force_velocity.py
)
set_tests_properties(thermal_force_velocity PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/thermal_force_velocity_work"
)

add_test(
    NAME thermal_force_slab
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_thermal_force_slab.py
)
set_tests_properties(thermal_force_slab PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/thermal_force_slab_work"
)

add_test(
    NAME thermal_force_tokamak
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_thermal_force_tokamak.py
)
set_tests_properties(thermal_force_tokamak PROPERTIES
    ENVIRONMENT
        "APPLETS_ROOT=${CMAKE_SOURCE_DIR};\
GORILLA_APPLETS_BIN=$<TARGET_FILE:gorilla_applets_main.x>;\
WORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/thermal_force_tokamak_work"
)

# Impurity tracing (i_option 17): charge-state sampling of marker_species_mod vs rate equations
add_executable(test_marker_species.x
    test_marker_species.f90
)
target_link_libraries(test_marker_species.x GORILLA_APPLETS)
target_include_directories(test_marker_species.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)
add_test(NAME marker_species COMMAND test_marker_species.x)

add_executable(test_marker_species_reject.x
    test_marker_species_reject.f90
)
target_link_libraries(test_marker_species_reject.x GORILLA_APPLETS)
target_include_directories(test_marker_species_reject.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)
add_test(NAME marker_species_reject COMMAND test_marker_species_reject.x)
set_tests_properties(marker_species_reject PROPERTIES
    PASS_REGULAR_EXPRESSION "ionisation of the bare nucleus")

# ADAS ADF11 reader (adf11_mod) and marker_species_mod adas_rates on a synthetic file
add_executable(test_adf11.x
    test_adf11.f90
)
target_link_libraries(test_adf11.x GORILLA_APPLETS)
target_include_directories(test_adf11.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)
add_test(NAME adf11 COMMAND test_adf11.x)

# Impurity tracing: straight-line neutral flight (neutral_flight_mod) on the analytic circular tokamak
add_executable(test_neutral_flight.x
    test_neutral_flight.f90
)
target_link_libraries(test_neutral_flight.x GORILLA_APPLETS)
target_include_directories(test_neutral_flight.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)
# case: n2, n_field_periods, sag tolerance for 'auto', expected 'auto' choice
foreach(nf_case "coarse;8;1;0.01;hybrid" "periods;6;3;0.01;hybrid" "fine;400;1;0.1;cartesian")
    list(GET nf_case 0 nf_name)
    list(SUBLIST nf_case 1 4 nf_args)
    set(nf_dir ${CMAKE_CURRENT_BINARY_DIR}/neutral_flight_${nf_name})
    file(MAKE_DIRECTORY ${nf_dir})
    foreach(f gorilla.inp field_divB0.inp preload_for_SYNCH.inp)
        configure_file(${CMAKE_SOURCE_DIR}/EXAMPLES/thermal_force_tokamak/${f} ${nf_dir}/${f} COPYONLY)
    endforeach()
    add_test(NAME neutral_flight_${nf_name} COMMAND test_neutral_flight.x ${nf_args} WORKING_DIRECTORY ${nf_dir})
endforeach()
# GORILLA stops with exit code 0 on input errors: pass only on the final message
set_tests_properties(neutral_flight_coarse neutral_flight_periods neutral_flight_fine PROPERTIES
    PASS_REGULAR_EXPRESSION "test_neutral_flight: all checks passed"
    FAIL_REGULAR_EXPRESSION "FAIL|ERROR")
