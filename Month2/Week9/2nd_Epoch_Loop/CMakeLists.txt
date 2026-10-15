cmake_minimum_required(VERSION 3.15)
project(SPPProcessor LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)
 
add_executable(spp_processor
    apps/spp_processor.cpp
    src/satellites.cpp
    src/gps_ls_solver.cpp
    src/utils.cpp
    src/rinex_obs.cpp)

# CORRECTED (2026-10-15): transmit-time version of the same app.
add_executable(spp_processor_crctd
    apps/spp_processor_crctd.cpp
    src/gps_ls_solver.cpp
    src/utils.cpp
    src/satellites.cpp
    src/rinex_obs.cpp)