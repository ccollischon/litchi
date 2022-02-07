cmake_minimum_required(VERSION 3.10)
# to use g++-10, if just g++ is too old: CC=gcc-10 CXX=g++-10 before calling cmake
# set the project name
project(Litchi)

if(DEFINED ENV{HEALPIX})
    set(healpix_parent_dir $ENV{HEALPIX} )
else()
    message(WARNING "HEALPIX not defined, should contain path to Healpix-folder (with include/, lib/ subdirectories).\n Set if not installed with package manager!")
endif()
set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)

#include_directories(${PROJECT_SOURCE_DIR}/pybind11)
add_subdirectory("pybind11")
include_directories(${healpix_parent_dir}/include/)
link_directories(${healpix_parent_dir}/lib/)

# add the executable
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=c++20 -fopenmp -fconcepts -O3")
set(CMAKE_CXX_STANDARD 20)


add_executable(litchi regnitz.cpp)
add_executable(tests tests.cpp)
target_link_libraries(litchi -lhealpix_cxx)
target_link_libraries(tests PUBLIC -lhealpix_cxx)

#find_package(pybind11 REQUIRED)
pybind11_add_module(litchieat pybindings.cpp)
target_link_libraries(litchieat PUBLIC -lhealpix_cxx)
