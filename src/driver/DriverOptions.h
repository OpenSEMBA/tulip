#pragma once

#include "SolverOptions.h"

#include <string>

namespace tulip {

struct DriverOptions {
	SolverOptions solverOptions;
	
	// Number of coefficients in the multipolar expansion 
	// for in-cell parameters calculations.
	int multipolarExpansionOrder{ 3 }; 

	bool exportParaViewSolution{ true };
	
	std::string exportFolder{ "./" };
};

}