cmake_minimum_required(VERSION 3.8)
project(rsa3d)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")       # -O3 in relase, in debug no optimisations
if(NOT APPLE)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp")  # use OpenMP
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp -static")    # use OpenMP
endif()
# openmp for clang on osx - see https://mac.r-project.org/openmp/ and download necessary libraries
if(APPLE)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Xclang -fopenmp")  # use OpenMP
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")    # use OpenMP
endif()
#if(UNIX AND NOT APPLE)
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static") # link dependencies statically (on Linux)
#endif()

add_subdirectory(statistics)    # statistics library
add_subdirectory(rsa3d)         # common rsa code and rsa executable
add_subdirectory(stat_test)     # statistical tests
add_subdirectory(unit_test)     # unit tests