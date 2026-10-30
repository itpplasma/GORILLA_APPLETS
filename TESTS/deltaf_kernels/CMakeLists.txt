add_executable(test_deltaf_kernels.x test_deltaf_kernels.f90)
target_link_libraries(test_deltaf_kernels.x
  GORILLA_APPLETS_DELTAF_KERNELS GORILLA)
target_include_directories(test_deltaf_kernels.x
  PUBLIC ${CMAKE_BINARY_DIR}/SRC_CORE)
add_test(NAME deltaf_kernels COMMAND test_deltaf_kernels.x
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
