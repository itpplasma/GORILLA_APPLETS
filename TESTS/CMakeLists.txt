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
