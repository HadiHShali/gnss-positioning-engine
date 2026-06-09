cmake_minimum_required(VERSION 3.15)
project(SatellitePosition LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
add_executable(satellite_position apps/satellite_position.cpp)
