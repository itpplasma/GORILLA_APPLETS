add_executable(test_multi_species_collision.x
    test_multi_species_collision.f90
)
target_link_libraries(test_multi_species_collision.x GORILLA_APPLETS)
target_include_directories(test_multi_species_collision.x
    PUBLIC ${CMAKE_BINARY_DIR}/SRC ${CMAKE_BINARY_DIR}/SRC_CORE)

add_test(NAME multi_species_collision COMMAND test_multi_species_collision.x)

add_executable(test_global_transport_fit.x
  test_global_transport_fit.f90
)

target_link_libraries(test_global_transport_fit.x
  GORILLA_APPLETS
)

target_include_directories(test_global_transport_fit.x PUBLIC
  ${CMAKE_BINARY_DIR}/SRC
  ${CMAKE_BINARY_DIR}/SRC_CORE
)

add_test(NAME test_global_transport_fit COMMAND test_global_transport_fit.x)

add_executable(test_kramers_moyal.x
  test_kramers_moyal.f90
)

target_link_libraries(test_kramers_moyal.x
  GORILLA_APPLETS
)

target_include_directories(test_kramers_moyal.x PUBLIC
  ${CMAKE_BINARY_DIR}/SRC
  ${CMAKE_BINARY_DIR}/SRC_CORE
)

add_test(NAME test_kramers_moyal COMMAND test_kramers_moyal.x)
