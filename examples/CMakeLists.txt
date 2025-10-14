# Pyramid flattening app
add_executable(openabf_example_basic_flattening src/BasicFlattening.cpp)
target_link_libraries(openabf_example_basic_flattening OpenABF::OpenABF)

add_executable(openabf_example_split_path src/SplitPath.cpp)
target_link_libraries(openabf_example_split_path OpenABF::OpenABF)

add_executable(openabf_example_flatten src/FlattenTool.cpp)
target_link_libraries(openabf_example_flatten OpenABF::OpenABF)

add_executable(openabf_example_clone src/CloneMesh.cpp)
target_link_libraries(openabf_example_clone OpenABF::OpenABF)