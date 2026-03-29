cmake_minimum_required(VERSION 3.12)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (GORILLA_APPLETS
  VERSION 1.0.0)

enable_language(Fortran)

if(UNIX AND NOT APPLE)
  include_directories(/usr/include)
elseif(APPLE)
  if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
    include_directories(/opt/local/include)
  elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    include_directories(/usr/local/include)
  endif()
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(netCDF REQUIRED)

enable_testing()

add_compile_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)
add_link_options(-g -fbacktrace -ffpe-trap=zero,overflow,invalid -fbounds-check -fopenmp)

if(APPLE)
  add_link_options(-L/opt/local/lib)
endif()

add_subdirectory(SRC)
add_subdirectory(SRC_CORE)
add_subdirectory(TESTS)

add_executable(gorilla_applets_main.x
  SRC/gorilla_applets_main.f90
)

target_link_libraries(gorilla_applets_main.x
  GORILLA_APPLETS
  )
target_include_directories(gorilla_applets_main.x PUBLIC ${CMAKE_BINARY_DIR}/SRC_CORE)
