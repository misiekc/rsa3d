cmake_minimum_required(VERSION 3.8)
project(rsa3d)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")   # -O3 in relase, in debug no optimisations
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp -static")    # link dependencies statically

add_subdirectory(statistics)
add_subdirectory(rsa3d)
add_subdirectory(stat_test)
add_subdirectory(unit_test)