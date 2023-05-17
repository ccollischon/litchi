cmake_minimum_required(VERSION 3.10)
# to use g++-10, if just g++ is too old: CC=gcc-10 CXX=g++-10 before calling cmake
# set the project name
project(Litchi)

if(DEFINED ENV{HEALPIX})
    set(healpix_parent_dir $ENV{HEALPIX} )
else()
    message(WARNING "Environment variable HEALPIX is not defined, should contain path to Healpix-folder (with include/, lib/ as subdirectories).\n Set it if Healpix was not installed with the package manager!")
endif()
set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)


add_subdirectory("pybind11")
include_directories(${healpix_parent_dir}/include/)
link_directories(${healpix_parent_dir}/lib/)

# add the executable
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=c++20 -fopenmp -fconcepts -O3")
set(CMAKE_CXX_STANDARD 20)

#litchieat as its own library, reducing redundant compilations
add_library(litchiwrap SHARED litchi_wrap.cpp)
target_link_libraries(litchiwrap PUBLIC -lhealpix_cxx)

#litchi itself
add_executable(litchi regnitz.cpp)
target_compile_definitions(litchi PUBLIC NDEBUG)
target_link_libraries(litchi litchiwrap -lhealpix_cxx)

#tests
add_executable(tests tests.cpp)
target_compile_definitions(tests PUBLIC THISRUNSINATEST)
target_link_libraries(tests PUBLIC litchiwrap -lhealpix_cxx)

add_executable(testGeometry testGeometry.cpp)
target_link_libraries(testGeometry PUBLIC litchiwrap -lhealpix_cxx)


#find_package(pybind11 REQUIRED)
pybind11_add_module(litchieat NO_EXTRAS pybindings.cpp)
target_compile_definitions(litchieat PUBLIC NDEBUG)
target_compile_definitions(litchieat PUBLIC THISISPYTHON)
target_link_libraries(litchieat PUBLIC litchiwrap -lhealpix_cxx)
