cmake_minimum_required(VERSION 3.16)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Building Release")
    set(CMAKE_BUILD_TYPE "Release")
else()
    message("Building Debug")
endif()

project(WiggleKiller)

set(WK_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(WK_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED)

add_subdirectory(src)