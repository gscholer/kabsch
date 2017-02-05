project(MyProjectUsingEigen)

cmake_minimum_required(VERSION 2.6.2)

find_package(Eigen3 3.2.0 REQUIRED )
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(EIGEN3_INCLUDE_DIR)

set(SRCS kabsch.cpp)
add_executable(kabsch ${SRCS})
