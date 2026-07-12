add_executable(test_multi_species_collision.x
    test_multi_species_collision.f90
)
target_link_libraries(test_multi_species_collision.x GORILLA_APPLETS)
target_include_directories(test_multi_species_collision.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)

add_test(NAME multi_species_collision COMMAND test_multi_species_collision.x)

add_executable(test_helical_core_profile.x
    test_helical_core_profile.f90
)
target_link_libraries(test_helical_core_profile.x GORILLA_APPLETS)
target_include_directories(test_helical_core_profile.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)
add_test(NAME helical_core_profile COMMAND test_helical_core_profile.x)

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
    NAME mono_energetic_transport
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_mono_energetic_transport.py
)
set_tests_properties(mono_energetic_transport PROPERTIES
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
