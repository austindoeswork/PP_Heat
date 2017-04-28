cmake_minimum_required(VERSION 3.6)
project(PP_Heat)

include_directories(/usr/include/mpi/)
find_package(MPI REQUIRED)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS -pthread)

set(SOURCE_FILES main.c)
add_executable(PP_Heat ${SOURCE_FILES})