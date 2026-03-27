#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

namespace pulmtln {

static const double VACUUM_RELATIVE_PERMITTIVITY{ 1.0 };

// Permittivity of Free Space (natural units)
static const double EPSILON0_NATURAL{ 1.0 };

// Permittivity of Free Space (units F/m)
static const double EPSILON0_SI{ 8.8541878176e-12 };

// Permeability of Free Space (natural units)
static const double MU0_NATURAL{ 1.0 };

// Permeability of Free Space (units H/m)
static const double MU0_SI{ 4.0e-7*M_PI};

}