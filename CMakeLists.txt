cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(DGHydro VERSION 0.1)

configure_file(DGConfig.h.in DGConfig.h)

# Add kokkos
add_subdirectory(src/kokkos build/kokkos)

add_executable(DGHydro)
add_subdirectory(src build)

target_include_directories(DGHydro PUBLIC
                                         "${PROJECT_BINARY_DIR}"
                                         )
target_include_directories(DGHydro PUBLIC
                                          src/kokkos/core/src
                                          src
                                         )

target_link_libraries(DGHydro Kokkos::kokkos)
