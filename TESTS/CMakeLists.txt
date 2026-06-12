add_executable(test_multi_species_collision.x
    test_multi_species_collision.f90
)
target_link_libraries(test_multi_species_collision.x GORILLA_APPLETS)
target_include_directories(test_multi_species_collision.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)

add_test(NAME multi_species_collision COMMAND test_multi_species_collision.x)

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
