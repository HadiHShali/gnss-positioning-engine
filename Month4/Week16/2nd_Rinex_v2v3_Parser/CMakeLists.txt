# ============================================================
# Week 16 Day 2 -- RINEX 2/3 version-aware parsers
# ============================================================
cmake_minimum_required(VERSION 3.15)
project(W16D2_Rinex_v2v3 LANGUAGES CXX)

# C++17 for everything in this folder (structured bindings, std::string features)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Headers live in include/  (rinex_obs.h, cycle_slip.h, satellites.h)
include_directories(include)

# MSVC only defines M_PI if this is set before <cmath> is included
# (check_station.cpp uses M_PI to convert the inclination to degrees)
add_definitions(-D_USE_MATH_DEFINES)

# ------------------------------------------------------------
# One static library with the shared code, compiled ONCE and
# linked into every app below:
#   rinex_obs.cpp   observation reader, RINEX 2 and 3
#   satellites.cpp  nav reader (RINEX 2 and 3), satellite positions,
#                   findBestEphemeris (skips unhealthy satellites)
# ------------------------------------------------------------
add_library(gnss_lib STATIC
    src/rinex_obs.cpp
    src/satellites.cpp
)

# ------------------------------------------------------------
# Apps
# ------------------------------------------------------------
# Week 15 cycle-slip scanner, now taking the obs file as an argument:
#   ./build/Release/test_phase_parser_allSat.exe data/<obs file>
add_executable(test_phase_parser_allSat apps/test_phase_parser_allSat.cpp)
target_link_libraries(test_phase_parser_allSat PRIVATE gnss_lib)

# One-station health check (obs + optional nav file):
#   ./build/Release/check_station.exe data/<obs file> [data/<nav file>]
add_executable(check_station apps/check_station.cpp)
target_link_libraries(check_station PRIVATE gnss_lib)