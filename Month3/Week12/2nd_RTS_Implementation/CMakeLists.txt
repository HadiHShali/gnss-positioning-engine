cmake_minimum_required(VERSION 3.15)
project(KFSPPProcessor LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)

# Library: KF class + reader + utils
add_library(kf_gps      STATIC src/kalman_filter_gps.cpp)
add_library(traj_reader STATIC src/trajectory_reader.cpp)
add_library(utils_lib   STATIC src/utils.cpp)

# Executable
add_executable(kf_spp_processor apps/kf_spp_processor.cpp)
target_link_libraries(kf_spp_processor PRIVATE
    kf_gps traj_reader utils_lib)

# CORRECTED (2026-10-24)
add_library(kf_gps_crctd STATIC src/kalman_filter_gps_crctd.cpp)
add_executable(kf_spp_processor_crctd apps/kf_spp_processor_crctd.cpp)
target_link_libraries(kf_spp_processor_crctd PRIVATE kf_gps_crctd traj_reader utils_lib)