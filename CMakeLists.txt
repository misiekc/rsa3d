cmake_minimum_required(VERSION 3.8)
project(rsa3d)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")       # -O3 in relase, in debug no optimisations
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")    # use OpenMP
if(UNIX AND NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static") # link dependencies statically (on Linux)
endif()

add_subdirectory(statistics)    # statistics library
add_subdirectory(rsa3d)         # common rsa code and rsa executable
add_subdirectory(stat_test)     # statistical tests
add_subdirectory(unit_test)     # unit tests