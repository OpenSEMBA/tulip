#pragma once

#include "FES.h"
#include "constants.h"

namespace pulmtln {

using multipolarCoefficient = std::pair<double, double>;
using multipolarCoefficients = std::vector<multipolarCoefficient>;

static double multipolarExpansion(
    const mfem::Vector& position,
    const multipolarCoefficients& ab,
    const mfem::Vector expansionCenter)
{
    // 2D multipolar expansion from:
    // TSOGTGEREL GANTUMUR, MULTIPOLE EXPANSIONS IN THE PLANE. 
    // 2016-04-16 lecture notes. 

    mfem::Vector rVec{ position };
    rVec -= expansionCenter;

    double r{ rVec.Norml2() };
    double phi{ std::atan2(rVec(1), rVec(0)) };

    double res{ 0.0 };
    for (int n{ 0 }; n < ab.size(); ++n) {
        const auto& an = ab[n].first;
        const auto& bn = ab[n].second;
        if (n == 0) {
            res -= an * std::log(r);
            assert(bn == 0.0); // b0 should always be zero.
        }
        else {
            res += (an * std::cos(n * phi) + bn * std::sin(n * phi)) / std::pow(r, n);
        }
    }
	res /= 2.0 * M_PI; 
    return res;
}

static double momentComponent(
    const mfem::Vector& position, int n, int component, const mfem::Vector& translation)
{
    assert(n > 0);
    assert(component == 0 || component == 1);

    if (n == 0) {
        return 1.0;
    }

    mfem::Vector pos{ position };
    pos -= translation;

    double r{ pos.Norml2() };
    double phi{ std::atan2(pos(1), pos(0)) };

    if (component == 0) {
        return std::pow(r, n) * std::cos(n * phi) / n;
    }
    else {
        return std::pow(r, n) * std::sin(n * phi) / n;
    }
}

}