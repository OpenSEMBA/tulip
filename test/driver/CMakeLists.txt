cmake_minimum_required(VERSION 3.21)

message(STATUS "Creating build system for driver_tests")

find_package(GTest CONFIG REQUIRED)

add_executable(driver_tests
    "DriverTest.cpp"
    "mfemTest.cpp"
    "ResultsTest.cpp"
    "DirectedGraphTest.cpp"
    "DomainTest.cpp"
    "multipolarExpansionTest.cpp"
    "ModelTest.cpp"
    "ParserTest.cpp"
    "SolverTest.cpp"
)

target_include_directories(driver_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(driver_tests
    tulip-driver
    GTest::gtest
    GTest::gtest_main
)

add_test(
    NAME driver_tests
    COMMAND driver_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

