#pragma once

#include <string>

namespace tulip {

struct SolverOptions {
	// Basis function order
	int order{3}; 
	bool printIterations{ false };
};

}