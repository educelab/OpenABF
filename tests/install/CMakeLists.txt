cmake_minimum_required(VERSION 3.15)
project(OpenABF_InstallTest)

find_package(OpenABF REQUIRED)

add_executable(install_test main.cpp)
target_link_libraries(install_test OpenABF::OpenABF)
