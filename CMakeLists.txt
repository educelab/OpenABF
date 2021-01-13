cmake_minimum_required(VERSION 3.17)
project(abf++)

set(CMAKE_CXX_STANDARD 14)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(abf_test main.cpp)