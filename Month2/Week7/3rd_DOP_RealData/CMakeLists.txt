cmake_minimum_required(VERSION 3.15)
project(Position_DOPSolver LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
include_directories("C:/eigen-3.4.0")
include_directories(include)
 
add_executable(position_dop_solver
    apps/real_data_position_withDOP.cpp
    src/satellites.cpp
    src/rinex_obs.cpp
    src/gps_ls_solver.cpp
    src/utils.cpp)
# CORRECTED (2026-10-15): transmit-time version of the same app.
add_executable(real_data_position_withDOP_crctd
    apps/real_data_position_withDOP_crctd.cpp
    src/gps_ls_solver.cpp
    src/utils.cpp
    src/satellites.cpp
    src/rinex_obs.cpp)