#pragma once

#include <string>

namespace tulip {

struct SolverOptions {
	// H1 basis function order used to solve the problem
	int order{3}; 

	// Prints iterations of the iterative solver.
	bool printIterations{ false };
};

}