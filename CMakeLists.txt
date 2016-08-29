cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(mcubes)

set(CMAKE_CXX_STANDARD 11)

include_directories("/usr/include/eigen3/")

add_definitions("-Wall")
set(SOURCES main.cpp vtk.cpp mesh.cpp)

add_executable(main ${SOURCES})

add_executable(cube cube.cpp)
