# CMakeLists.txt for NMEA parser version 2
cmake_minimum_required(VERSION 3.15)
project(NMEAParserV2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(nmea_parser_v2 apps/nmea_parser_v2.cpp)