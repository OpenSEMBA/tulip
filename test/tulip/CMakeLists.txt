cmake_minimum_required(VERSION 3.21)

message(STATUS "Creating build system for ${PROJECT_NAME}_tests")

find_package(GTest CONFIG REQUIRED)

add_executable(${PROJECT_NAME}_tests
    "ElectrostaticSolverTest.cpp"
    "mfemTest.cpp"
    "DriverTest.cpp"
    "ParserTest.cpp"
    "ResultsTest.cpp"
    "ModelTest.cpp"
    "DirectedGraphTest.cpp"
    "DomainTest.cpp"
    "multipolarExpansionTest.cpp"
)

target_include_directories(${PROJECT_NAME}_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME}_tests
    pulmtln-core
    GTest::gtest
    GTest::gtest_main
)

