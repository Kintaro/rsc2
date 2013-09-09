cmake_minimum_required(VERSION 2.6)
project(sash CXX)

# set mpi compiler
include_directories(${Boost_INCLUDE_DIRS})
include_directories(../../)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS})

add_definitions(-std=c++0x -g -Wall -fPIC -Wno-attributes)
add_library(sash SHARED Sash.cpp ../../FileUtil.cpp)

# link with mpi
target_link_libraries(sash ${Boost_LIBRARIES})

set_target_properties(sash PROPERTIES PREFIX "lib")
set_target_properties(sash PROPERTIES SUFFIX ".extension")
