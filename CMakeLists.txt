cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(DGHydro VERSION 0.1)
add_executable(DGHydro src/main.cpp)

configure_file(DGConfig.h.in DGConfig.h)
target_include_directories(DGHydro PUBLIC
                                         "${PROJECT_BINARY_DIR}"
                                         )
