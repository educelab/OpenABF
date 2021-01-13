cmake_minimum_required(VERSION 3.17)
project(OpenABF)

set(CMAKE_CXX_STANDARD 14)

include_directories(${CMAKE_SOURCE_DIR}/include)
set(public_hdrs
    include/OpenABF/ABF.hpp
    include/OpenABF/HalfEdgeMesh.hpp
)

add_executable(abf_test ${public_hdrs} main.cpp)