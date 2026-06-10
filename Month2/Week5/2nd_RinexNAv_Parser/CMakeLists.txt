# CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(RinexNavParser LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
add_executable(rinex_nav_parser apps/rinex_nav_parser.cpp)
