cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(b3d_fromcascade CXX)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_b3d COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/b3d/software)
add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/b3d/software/include")

include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/b3d/software/include ${GITHOME_MSU}/commonutils/software/include)
link_directories(${GITHOME_MSU}/b3d/software/lib ${GITHOME_MSU}/commonutils/software/lib)
message("Libraries: b3d ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

find_library(${GITHOME_MSU}/b3d/software/lib b3d)
find_library(${GITHOME_MSU}/commonutils/software/lib commonutils)

set(SOURCE_FILES ./b3dmain.cc)
add_executable(b3dtest ${SOURCE_FILES})
add_dependencies(b3dtest extern_b3d extern_commonutils)
target_link_libraries(b3dtest b3d commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY} ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GITHOME_MSU}/b3d/software/lib/libb3d.a)
