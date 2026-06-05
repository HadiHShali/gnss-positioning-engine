# CMakeLists.txt for RINEX Header Parser
cmake_minimum_required(VERSION 3.15)
project(RinexHeaderParser LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
add_executable(rinex_header_parser apps/rinex_header_parser.cpp)
