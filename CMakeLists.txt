cmake_minimum_required(VERSION 3.10)

# set the project name
project(Litchi)

set(healpix_parent_dir "/userdata/data/collischon/Healpix_3.70/")
include_directories(${healpix_parent_dir}include/)
link_directories(${healpix_parent_dir}lib/)

# add the executable
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=c++20 -fopenmp -fconcepts -O3")


set(CMAKE_CXX_STANDARD 20)
add_executable(regnitz regnitz.cpp)
target_link_libraries(regnitz -lhealpix_cxx)
