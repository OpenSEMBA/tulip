cmake_minimum_required(VERSION 3.21)
project(step2gmsh)

set(CMAKE_CXX_STANDARD 17)

# -----------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------
find_package(gmsh CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

# -----------------------------------------------------------------------
# Library
# -----------------------------------------------------------------------
set(STEP2GMSH_LIB step2gmsh-core)

add_library(${STEP2GMSH_LIB} STATIC
    src/Graph.cpp
    src/BoundingBox.cpp
    src/ShapesClassification.cpp
    src/AreaExporterService.cpp
    src/Mesher.cpp
)

target_include_directories(${STEP2GMSH_LIB} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(${STEP2GMSH_LIB} PUBLIC
    $<IF:$<TARGET_EXISTS:gmsh::shared>,gmsh::shared,gmsh::lib>
    nlohmann_json::nlohmann_json
)

# -----------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------
enable_testing()

find_package(GTest CONFIG REQUIRED)

set(STEP2GMSH_TESTS step2gmsh_tests)

add_executable(${STEP2GMSH_TESTS}
    test/GraphTest.cpp
    test/BoundingBoxTest.cpp
    test/ShapesClassificationTest.cpp
    test/MesherTest.cpp
)

target_include_directories(${STEP2GMSH_TESTS} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/test
)

target_link_libraries(${STEP2GMSH_TESTS}
    ${STEP2GMSH_LIB}
    GTest::gtest
    GTest::gtest_main
)

# Run tests from the project root so ./testData/ resolves correctly
add_test(
    NAME ${STEP2GMSH_TESTS}
    COMMAND ${STEP2GMSH_TESTS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/..
)
