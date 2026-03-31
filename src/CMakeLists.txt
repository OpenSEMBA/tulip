cmake_minimum_required(VERSION 3.21)

set(EXECUTABLE_NAME "tulip")
message(STATUS "Creating build system for ${EXECUTABLE_NAME}")

add_executable(${EXECUTABLE_NAME}
    "tulip.cpp"
)

target_include_directories(
    "adapter/",
    "driver/"
)

target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${LIBRARY_NAME} Boost::program_options)