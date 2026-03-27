cmake_minimum_required(VERSION 3.21)

set(LIBRARY_NAME "tulip-driver")
message(STATUS "Creating build system for ${LIBRARY_NAME}")

find_package(nlohmann_json CONFIG REQUIRED)

find_package(Boost REQUIRED graph program_options)

add_library(${LIBRARY_NAME} STATIC
    "FES.cpp"
    "Solver.cpp"
    "Driver.cpp"
    "Parser.cpp"
    "Model.cpp"
    "Materials.cpp"
    "Results.cpp"
    "DirectedGraph.cpp"
    "Domain.cpp"
)

target_include_directories(${LIBRARY_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(${LIBRARY_NAME} PUBLIC mfem nlohmann_json::nlohmann_json Boost::graph)
