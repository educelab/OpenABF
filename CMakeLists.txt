cmake_minimum_required(VERSION 3.15...3.30 FATAL_ERROR)
project(OpenABF VERSION 2.0.0)

# Modules
include(FetchContent)
include(CMakeDependentOption)
include(GNUInstallDirs)

# Setup project directories
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Find dependencies
find_package(Eigen3 5 QUIET CONFIG)
if(NOT Eigen3_FOUND)
    find_package(Eigen3 3.3 REQUIRED CONFIG)
endif()

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
        include/OpenABF/MeshIO.hpp
        include/OpenABF/MeshIOUtils.hpp
        include/OpenABF/MeshIOFormats.hpp
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
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(OpenABF INTERFACE cxx_std_17)
set_target_properties(OpenABF
    PROPERTIES
        PUBLIC_HEADER "${public_hdrs}"
)
target_link_libraries(OpenABF
    INTERFACE
        Eigen3::Eigen
)

# Updated all destinations to use standard directory variables
install(
    TARGETS OpenABF
    EXPORT OpenABFTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenABF
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/OpenABF
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