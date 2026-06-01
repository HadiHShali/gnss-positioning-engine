# CMakeLists.txt for NMEA Parser
cmake_minimum_required(VERSION 3.15)
project(NMEAParser LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
# Add the executable
add_executable(nmea_parser apps/nmea_parser.cpp)
