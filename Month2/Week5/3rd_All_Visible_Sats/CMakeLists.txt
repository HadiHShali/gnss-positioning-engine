cmake_minimum_required(VERSION 3.15)
project(AllVisibleSats LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
add_executable(all_visible_sats apps/all_visible_sats.cpp)
