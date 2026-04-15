cmake_minimum_required(VERSION 3.21)

find_package(Boost REQUIRED COMPONENTS program_options)

set(EXECUTABLE_NAME "tulip")
message(STATUS "Creating build system for ${EXECUTABLE_NAME}")

add_executable(${EXECUTABLE_NAME}
    "tulip.cpp"
    "Launcher.cpp"
)

target_include_directories(${EXECUTABLE_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/adapter
    ${CMAKE_CURRENT_SOURCE_DIR}/driver
)

target_link_libraries(${EXECUTABLE_NAME} PRIVATE 
    tulip-adapter
    tulip-driver
    Boost::program_options
)