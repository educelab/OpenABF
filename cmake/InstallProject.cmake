include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    OpenABFConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion
    ARCH_INDEPENDENT
)
configure_file(cmake/OpenABFConfig.cmake.in OpenABFConfig.cmake @ONLY)

install(
    EXPORT OpenABFTargets
    FILE OpenABFTargets.cmake
    NAMESPACE OpenABF::
    DESTINATION lib/cmake/OpenABF
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/OpenABFConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/OpenABFConfigVersion.cmake"
    DESTINATION lib/cmake/OpenABF
)