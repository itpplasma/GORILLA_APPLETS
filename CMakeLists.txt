cmake_minimum_required(VERSION 3.12)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (GORILLA_APPLETS
  VERSION 1.0.0)

enable_language(Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(netCDF REQUIRED)
find_package(netCDF-Fortran REQUIRED)

get_target_property(_NETCDFF_INC netCDF::netcdff INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${netCDF_INCLUDE_DIR} ${_NETCDFF_INC})

add_compile_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)
add_link_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)

add_subdirectory(SRC)
add_subdirectory(SRC_CORE)

add_executable(gorilla_applets_main.x
  SRC/gorilla_applets_main.f90
)

target_link_libraries(gorilla_applets_main.x
  GORILLA_APPLETS
  )
target_include_directories(gorilla_applets_main.x PUBLIC ${CMAKE_BINARY_DIR}/SRC_CORE)