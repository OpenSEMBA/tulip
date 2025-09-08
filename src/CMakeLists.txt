cmake_minimum_required(VERSION 3.21)

set(LIBRARY_NAME "pulmtln-core")
message(STATUS "Creating build system for ${LIBRARY_NAME}")

find_package(nlohmann_json CONFIG REQUIRED)

find_package(Boost REQUIRED graph program_options)
include_directories(${Boost_INCLUDE_DIRS})

add_library(${LIBRARY_NAME} STATIC 
	"FES.cpp"
	"ElectrostaticSolver.cpp"
	"Driver.cpp" 
	"Parser.cpp"
	"Model.cpp"
	"Materials.cpp"
	"Results.cpp" 
	"DirectedGraph.cpp"
	"Domain.cpp")

target_link_libraries(${LIBRARY_NAME} mfem nlohmann_json Boost::graph)

set(EXECUTABLE_NAME "pulmtln")
message(STATUS "Creating build system for ${EXECUTABLE_NAME}")
add_executable(${EXECUTABLE_NAME} 
	"pulmtln.cpp"
)


target_link_libraries(${EXECUTABLE_NAME} ${LIBRARY_NAME} Boost::program_options)