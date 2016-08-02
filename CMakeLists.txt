cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(mcubes)

set(CMAKE_CXX_STANDARD 11)

set(SOURCES main.cpp vtk.cpp)

add_executable(main ${SOURCES})
