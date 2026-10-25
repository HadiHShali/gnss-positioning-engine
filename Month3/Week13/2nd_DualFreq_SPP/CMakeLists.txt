cmake_minimum_required(VERSION 3.15)
project(SPPProcessorDualFreq LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)

# Library: satellite ephemeris/atmospheric models + LS solver + utils + RINEX parser
add_library(satellites_lib STATIC src/satellites.cpp)
add_library(solver_lib     STATIC src/gps_ls_solver.cpp)
add_library(utils_lib      STATIC src/utils.cpp)
add_library(rinex_obs_lib  STATIC src/rinex_obs.cpp)

# Executable
add_executable(spp_processor_dualFreq apps/spp_processor_dualFreq.cpp)
target_link_libraries(spp_processor_dualFreq PRIVATE
    satellites_lib solver_lib utils_lib rinex_obs_lib)

# CORRECTED (2026-10-24): transmit-time version of the same app.
add_executable(spp_processor_dualFreq_crctd
    apps/spp_processor_dualFreq_crctd.cpp
    src/gps_ls_solver.cpp
    src/utils.cpp
    src/satellites.cpp
    src/rinex_obs.cpp)