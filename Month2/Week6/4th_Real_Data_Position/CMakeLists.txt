cmake_minimum_required(VERSION 3.15)
project(RealDataPosition LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
include_directories("C:/eigen-3.4.0")
include_directories(include)
 
add_executable(real_data_position
    apps/real_data_position.cpp
    src/gps_ls_solver.cpp
    src/utils.cpp
    src/satellites.cpp
    src/rinex_obs.cpp)

# CORRECTED (2026-10-12): transmit-time version of the same app.
# Same library sources; only the app file differs.
add_executable(real_data_position_crctd
    apps/real_data_position_crctd.cpp
    src/gps_ls_solver.cpp
    src/utils.cpp
    src/satellites.cpp
    src/rinex_obs.cpp)