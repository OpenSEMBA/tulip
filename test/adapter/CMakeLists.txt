cmake_minimum_required(VERSION 3.21)

message(STATUS "Creating build system for adapter_tests")

find_package(GTest CONFIG REQUIRED)

add_executable(adapter_tests
    GraphTest.cpp
    BoundingBoxTest.cpp
    ShapesClassificationTest.cpp
    MesherTest.cpp
)

target_include_directories(adapter_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(adapter_tests
    tulip-adapter
    GTest::gtest
    GTest::gtest_main
)

add_test(
    NAME adapter_tests
    COMMAND adapter_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
