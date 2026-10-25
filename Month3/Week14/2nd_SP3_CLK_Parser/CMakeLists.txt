cmake_minimum_required(VERSION 3.15)
project(SP3ClkParser LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
include_directories("C:/eigen-3.4.0")
include_directories(include)
add_definitions(-D_USE_MATH_DEFINES)
 
# Library: SP3 + CLK parsers (interpolation.h is header-only, no .cpp needed)
add_library(sp3_lib STATIC src/sp3_reader.cpp)
add_library(clk_lib STATIC src/clk_reader.cpp)
 
# Executable: standalone test program
add_executable(test_sp3_clk apps/test_sp3_clk.cpp)
target_link_libraries(test_sp3_clk PRIVATE sp3_lib clk_lib)

# CORRECTED (2026-10-24): same test app, corrected CLK reader
add_executable(test_sp3_clk_crctd apps/test_sp3_clk.cpp
    src/sp3_reader.cpp src/clk_reader_crctd.cpp)