#pragma once

#include <string>

namespace pulmtln {

struct SolverOptions {
	// Basis function order
	int order{3}; 
	bool printIterations{ false };
};

}