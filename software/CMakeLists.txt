cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})
project(b3d CXX C)
message("CXX is ${CXX}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++17 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include)

file(GLOB b3d_src src/*.cc)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(b3d STATIC ${b3d_src})
target_link_libraries(b3d ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Boost_LIBRARIES})
message("Libraries: ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}")




