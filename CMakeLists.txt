cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg --coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -march=native")

project(fastbe VERSION 1.0)

add_subdirectory(src/)
add_subdirectory(tests/)

target_include_directories(fastbe PUBLIC "${fastbe_SOURCE_DIR}/include")
target_include_directories(fastbe PUBLIC "${fastbe_SOURCE_DIR}/third-party/csv-parser/include")

# turn off JSON library tests
set(JSON_BuildTests OFF CACHE INTERNAL "")

# compile libraries
add_subdirectory(third-party/spdlog)
add_subdirectory(third-party/pprint)
add_subdirectory(third-party/csv-parser)
add_subdirectory(third-party/argparse)
add_subdirectory(third-party/json)
add_subdirectory(third-party/gtest)
add_subdirectory(third-party/nanobench)

# add include directories for all tests
target_include_directories(piecewiselinear_test PUBLIC "${fastbe_SOURCE_DIR}/include")
target_include_directories(piecewiselinear_test PUBLIC "${fastbe_SOURCE_DIR}/third-party/csv-parser/include")

enable_testing()
