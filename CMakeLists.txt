cmake_minimum_required(VERSION 3.16)

project(HOURGLASS_PLUS)

option(ADAPT "Enable adaptivity" ON)
option(BUILD_EXAMPLES "Build the examples" ON)
option(BUILD_BENCHMARKS "Build the benchmark targets" ON)

set(CMAKE_CXX_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Ofast")
endif ()

if (ADAPT)
  add_definitions(-DADAPT)
endif ()

find_package(OpenSSL REQUIRED)

add_library(hourglasslib
    hashutil.cc
)
target_include_directories(hourglasslib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(hourglasslib PUBLIC OpenSSL::SSL OpenSSL::Crypto)

if (BUILD_EXAMPLES)
    message(STATUS "Building examples")
    add_executable(hourglass example.cpp)
    target_link_libraries(hourglass hourglasslib)
    add_executable(hourglass_plus example_plus.cpp)
    target_link_libraries(hourglass_plus hourglasslib)
endif ()

if (BUILD_BENCHMARKS)
    message(STATUS "Building benchmarks")
    add_subdirectory(bench)
endif ()
