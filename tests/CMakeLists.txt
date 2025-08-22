## Google Test ##
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.17.0
        EXCLUDE_FROM_ALL
)
set(INSTALL_GTEST OFF CACHE BOOL OFF FORCE)
FetchContent_MakeAvailable(googletest)

## Build the tests ##
set(tests
    src/TestHalfEdgeMesh.cpp
    src/TestMath.cpp
    src/TestVec.cpp
)

foreach(src ${tests})
    get_filename_component(filename ${src} NAME_WE)
    set(testname OpenABF_${filename})
    add_executable(${testname} ${src})
    target_link_libraries(${testname}
        OpenABF::OpenABF
        gtest_main
        gmock_main
    )
    add_test(
        NAME ${testname}
        WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
        COMMAND ${testname}
    )
endforeach()