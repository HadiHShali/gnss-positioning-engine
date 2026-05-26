# CMakeLists.txt: Describes how to build the KF library project

# SECTION1: Project Setup ------------------------------------
cmake_minimum_required(VERSION 3.15)
project(2nd_KF_Library_Cmake LANGUAGES CXX)

# SECTION2: C++ Standard -------------------------------------
set(CMAKE_CXX_Standard 17)
Set(CMAKE_CXX_STANDARD_REQUIRED ON)

# SECTION3: Tell Compiler where to find header files ---------
# CMAKE_SOURCE_DIR is auotomatically set to the folder of this CMakeLists.txt
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(C:/eigen-3.4.0)

# SECTION4: Define the executables ---------------------------
# This creates app_kf_1d.exe from these two .cpp files
add_executable(app_kf_1d
	src/KalmanFilter.cpp
	apps/app_kf_1d.cpp
)

# SECTION5: Helpful message after configuring -----------------
message(STATUS "KF Library configured successfully")