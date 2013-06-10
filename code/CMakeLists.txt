cmake_minimum_required(VERSION 2.6)
# set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# set mpi compiler
include(CMakeForceCompiler)
CMAKE_FORCE_CXX_COMPILER(mpic++ "MPI C++ Compiler")

project(rsc2 CXX)

find_package(Boost 1.4 COMPONENTS mpi serialization REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(src)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS})

add_definitions(-std=c++0x -Wall -Wextra -Werror -pedantic -Wno-unused-parameter)
add_executable(rsc2 src/Daemon.cpp src/Options.cpp src/DefaultOptions.cpp src/FileUtil.cpp src/VecDataBlock.cpp src/AbstractSetManager.cpp src/rsc2.cpp src/RscClusterer.cpp)
# link with mpi
target_link_libraries(rsc2 ${Boost_LIBRARIES})
add_subdirectory(src/index_structures/sash)
