#pragma once

#include "SolverOptions.h"

#include <string>

namespace tulip {

struct DriverOptions {
	SolverOptions solverOptions;
	
	// Number of coefficients in the multipolar expansion 
	// for in-cell parameters calculations.
	int multipolarExpansionOrder{ 3 }; 

	// Exports Paraview solution.
	bool exportParaViewSolution{ true };
	
	// Sets export folder for all files relative to the input file.
	std::string exportFolder{ "./" };
};

}