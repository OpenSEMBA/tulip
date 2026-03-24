#pragma once

#include <string>

namespace pulmtln {

struct SolverOptions {
	int order{3};     // Basis function order
	bool printIterations{ false };
};

struct DriverOptions {
	SolverOptions solverOptions;
	
	// Number of coefficients in the multipolar expansion 
	// for in-cell parameters calculations.
	int multipolarExpansionOrder{ 5 }; 

	bool exportParaViewSolution{ true };
	
	std::string exportFolder{ "./" };
};

}