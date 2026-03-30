add_executable(test_multi_species_collision.x
    test_multi_species_collision.f90
)
target_link_libraries(test_multi_species_collision.x GORILLA_APPLETS)
target_include_directories(test_multi_species_collision.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)

add_test(NAME multi_species_collision COMMAND test_multi_species_collision.x)
