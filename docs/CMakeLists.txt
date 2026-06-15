# Set Doxygen options
set(DOXYGEN_DOT_IMAGE_FORMAT svg)
set(DOXYGEN_EXTRACT_STATIC YES)
set(DOXYGEN_EXTRACT_PRIVATE YES)
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_GENERATE_LATEX NO)
set(DOXYGEN_HTML_EXTRA_STYLESHEET ${PROJECT_SOURCE_DIR}/docs/style/doxygen_extra.css)
set(DOXYGEN_SHOW_FILES NO)
set(DOXYGEN_SHOW_USED_FILES NO)
set(DOXYGEN_SOURCE_BROWSER NO)
set(DOXYGEN_VERBATIM_HEADERS NO)
set(DOXYGEN_TEMPLATE_RELATIONS YES)
set(DOXYGEN_QUIET YES)

set(DOXYGEN_CITE_BIB_FILES ${PROJECT_SOURCE_DIR}/docs/citations.bib)
set(DOXYGEN_IMAGE_PATH ${PROJECT_SOURCE_DIR}/docs/images/)

# Build list of example applications
set(DOXYGEN_EXAMPLE_PATH ${PROJECT_SOURCE_DIR}/examples/src/)
set(EXAMPLE_FILES
    ${PROJECT_SOURCE_DIR}/examples/src/BasicFlattening.cpp
    ${PROJECT_SOURCE_DIR}/examples/src/SplitPath.cpp
    ${PROJECT_SOURCE_DIR}/examples/src/MultiChartFlatten.cpp
    ${PROJECT_SOURCE_DIR}/examples/src/FlattenTool.cpp
    ${PROJECT_SOURCE_DIR}/examples/src/CloneMesh.cpp
    ${PROJECT_SOURCE_DIR}/examples/src/AlternativeSolvers.cpp
    ${PROJECT_SOURCE_DIR}/examples/src/BenchmarkFlattening.cpp
)

set(DOXYGEN_PROJECT_NAME "OpenABF")
set(DOXYGEN_PROJECT_LOGO ${PROJECT_SOURCE_DIR}/docs/images/logo.svg)
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${PROJECT_SOURCE_DIR}/docs/index.md")

set(DOXYGEN_STRIP_FROM_INC_PATH ${OPENABF_LINK_DIR})
set(DOXYGEN_STRIP_FROM_PATH ${OPENABF_LINK_DIR})

doxygen_add_docs(
    docs
        ${PROJECT_SOURCE_DIR}/docs/index.md
        ${PROJECT_SOURCE_DIR}/docs/installation.md
        ${PROJECT_SOURCE_DIR}/docs/modules.dox
        ${OPENABF_LINK_DIR}
        ${EXAMPLE_FILES}
        COMMENT "Generate Doxygen documentation"
)