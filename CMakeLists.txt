#
# This file is part of litchi, a lightweight C++ library
# for Minkowski analysis
#
# Copyright (C) 2021-2025 Caroline Collischon <caroline.collischon@fau.de>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
#


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

#include_directories(${PROJECT_SOURCE_DIR}/pybind11)
add_subdirectory("pybind11")
include_directories(${healpix_parent_dir}/include/)
link_directories(${healpix_parent_dir}/lib/)

# add the executable
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++20 -fopenmp -fconcepts")
set(CMAKE_CXX_STANDARD 20)


add_executable(litchi regnitz.cpp)
target_compile_options(litchi PUBLIC "-O3")
target_compile_definitions(litchi PUBLIC NDEBUG)
target_link_libraries(litchi -lhealpix_cxx)


add_executable(litchi_debug EXCLUDE_FROM_ALL regnitz.cpp)
target_compile_options(litchi_debug PUBLIC "-pg")
target_compile_options(litchi_debug PUBLIC "-g")
target_compile_options(litchi_debug PUBLIC "-O2")
target_link_libraries(litchi_debug -lhealpix_cxx)
target_link_options(litchi_debug PUBLIC "-pg")
target_link_options(litchi_debug PUBLIC "-g")

add_executable(tests EXCLUDE_FROM_ALL tests.cpp)
target_compile_options(tests PUBLIC "-pg")
target_compile_options(tests PUBLIC "-g")
target_compile_options(tests PUBLIC "-O2")
target_compile_definitions(tests PUBLIC THISRUNSINATEST)
target_link_libraries(tests PUBLIC -lhealpix_cxx)
target_link_options(tests PUBLIC "-pg")
target_link_options(tests PUBLIC "-g")

add_executable(testGeometry EXCLUDE_FROM_ALL testGeometry.cpp)
target_link_libraries(testGeometry PUBLIC -lhealpix_cxx)


#find_package(pybind11 REQUIRED)
pybind11_add_module(litchieat NO_EXTRAS pybindings.cpp)
target_compile_options(litchieat PUBLIC "-O3")
target_compile_definitions(litchieat PUBLIC NDEBUG)
target_compile_definitions(litchieat PUBLIC THISISPYTHON)
target_link_libraries(litchieat PUBLIC -lhealpix_cxx)


pybind11_add_module(litchieat_debug EXCLUDE_FROM_ALL NO_EXTRAS pybindings_debug.cpp)
target_compile_definitions(litchieat_debug PUBLIC THISISPYTHON)
target_link_libraries(litchieat_debug PUBLIC -lhealpix_cxx)
target_compile_options(litchieat_debug PUBLIC "-pg")
target_compile_options(litchieat_debug PUBLIC "-g")
