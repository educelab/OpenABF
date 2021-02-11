cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(OpenABF VERSION 1.0)

# Setup project directories
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Modules
include(FetchContent)
include(CMakeDependentOption)

# Find dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

option(OPENABF_MULTIHEADER "Build the multi-header version of OpenABF" off)

# OpenABF Library
if(OPENABF_MULTIHEADER)
    set(public_hdrs
        include/OpenABF/OpenABF.hpp
        include/OpenABF/ABF.hpp
        include/OpenABF/ABFPlusPlus.hpp
        include/OpenABF/Exceptions.hpp
        include/OpenABF/HalfEdgeMesh.hpp
        include/OpenABF/AngleBasedLSCM.hpp
        include/OpenABF/Math.hpp
        include/OpenABF/Vec.hpp
    )
    set(OPENABF_LINK_DIR ${PROJECT_SOURCE_DIR}/include/)
    message(STATUS "Using multi-header OpenABF library: ${OPENABF_LINK_DIR}")
else()
    set(public_hdrs single_include/OpenABF/OpenABF.hpp)
    set(OPENABF_LINK_DIR ${PROJECT_SOURCE_DIR}/single_include/)
    message(STATUS "Using single-header OpenABF library: ${OPENABF_LINK_DIR}")
endif()

add_library(OpenABF INTERFACE)
add_library(OpenABF::OpenABF ALIAS OpenABF)
target_include_directories(OpenABF
    INTERFACE
        $<BUILD_INTERFACE:${OPENABF_LINK_DIR}>
        $<INSTALL_INTERFACE:include>
)
target_compile_features(OpenABF INTERFACE cxx_std_14)
set_target_properties(OpenABF
    PROPERTIES
        PUBLIC_HEADER "${public_hdrs}"
)
target_link_libraries(OpenABF
    INTERFACE
        Eigen3::Eigen
)
install(
    TARGETS OpenABF
    EXPORT OpenABFTargets
    ARCHIVE DESTINATION "lib"
    LIBRARY DESTINATION "lib"
    INCLUDES DESTINATION "include/OpenABF"
    PUBLIC_HEADER DESTINATION "include/OpenABF"
)

# Docs
find_package(Doxygen OPTIONAL_COMPONENTS dot)
CMAKE_DEPENDENT_OPTION(OPENABF_BUILD_DOCS "Build Doxygen documentation" off "DOXYGEN_FOUND" off)
if(OPENABF_BUILD_DOCS)
    add_subdirectory(docs)
endif()

# Examples
option(OPENABF_BUILD_EXAMPLES "Compile OpenABF example applications" off)
if(OPENABF_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Tests
option(OPENABF_BUILD_TESTS "Compile OpenABF unit tests" off)
if(OPENABF_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Install
include(InstallProject)