# ============================================================================
#  CMakeLists.txt  --  Month 4, Week 15, Day 2 : phase parser + cycle slips
# ----------------------------------------------------------------------------
#  Folder layout this file expects:
#      include/   rinex_obs.h, cycle_slip.h
#      src/       rinex_obs.cpp
#      apps/      test_phase_parser.cpp
#  Build, from this folder:
#      cmake -S . -B build
#      cmake --build build
# ============================================================================


# ============================================================================
#  SECTION 1 -- PROJECT AND LANGUAGE STANDARD
# ----------------------------------------------------------------------------
#  cmake_minimum_required : the oldest CMake version allowed to run this file.
#  project                : the project's name, and its language (C++).
#  CMAKE_CXX_STANDARD 17  : compile as C++17.
#  ..._REQUIRED ON        : stop with an error if the compiler cannot do
#                           C++17, instead of quietly using an older standard.
# ============================================================================
cmake_minimum_required(VERSION 3.15)
project(PhaseParser LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# ============================================================================
#  SECTION 2 -- WHERE THE HEADERS ARE, AND COMPILER SETTINGS
# ----------------------------------------------------------------------------
#  include_directories(include) : lets every .cpp write
#                           #include "rinex_obs.h"  without a folder name.
#  Eigen path             : not used by today's files. Kept for later days.
#  _USE_MATH_DEFINES      : makes M_PI available with the Microsoft compiler.
#                           Not used today; harmless.
#  Warnings               : ask the compiler to report suspicious code.
#                           /W4 is the Microsoft spelling, -Wall -Wextra the
#                           GCC / Clang spelling. The if() picks the right one.
# ============================================================================
include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra)
endif()


# ============================================================================
#  SECTION 3 -- THE LIBRARY : THE RINEX READER
# ----------------------------------------------------------------------------
#  add_library(... STATIC ...) : compile rinex_obs.cpp once, into a library.
#                           Every program that needs the reader links to it,
#                           so the reader is never compiled twice.
#  cycle_slip.h is not listed anywhere. It is header-only (constants and
#  inline functions), so there is nothing to compile for it.
# ============================================================================
add_library(rinex_obs_lib STATIC src/rinex_obs.cpp)


# ============================================================================
#  SECTION 4 -- THE PROGRAM
# ----------------------------------------------------------------------------
#  add_executable        : build test_phase_parser from its .cpp file.
#  target_link_libraries : attach the reader library to it. PRIVATE means
#                          "this program uses the library", nothing more.
# ============================================================================
add_executable(test_phase_parser apps/test_phase_parser.cpp)
target_link_libraries(test_phase_parser PRIVATE rinex_obs_lib)