cmake_minimum_required(VERSION 3.21)

set(ADAPTER_LIB "tulip-adapter")
message(STATUS "Creating build system for ${ADAPTER_LIB}")

find_package(gmsh CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

add_library(${ADAPTER_LIB} STATIC
    Graph.cpp
    BoundingBox.cpp
    ShapesClassification.cpp
    AreaExporterService.cpp
    Adapter.cpp
)

target_include_directories(${ADAPTER_LIB} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${ADAPTER_LIB} PUBLIC
    $<IF:$<TARGET_EXISTS:gmsh::shared>,gmsh::shared,gmsh::lib>
    nlohmann_json::nlohmann_json
)
