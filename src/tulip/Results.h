#pragma once

#include <nlohmann/json.hpp>

#include "FES.h"
#include "Domain.h"
#include "multipolarExpansion.h"

namespace pulmtln {


void saveToJSONFile(const nlohmann::json&, const std::string& filename);

struct PULParameters {
    PULParameters() = default;
    PULParameters(const nlohmann::json&);

    bool operator==(const PULParameters&) const;

    mfem::DenseMatrix getCapacitiveCouplingCoefficients() const;

    nlohmann::json toJSON() const;

    mfem::DenseMatrix L, C; // Stored in SI units.
};

struct PULParametersByDomain {
    DomainTree domainTree;
    std::map<Domain::Id, PULParameters> domainToPUL;
};

struct FloatingPotentials {
    mfem::DenseMatrix electric, magnetic;
};

struct FieldReconstruction {
    double innerRegionAveragePotential;
    std::array<double,2> expansionCenter;
    multipolarCoefficients ab; // Stored in natural units.
    std::map<MaterialId, double> conductorPotentials;

    nlohmann::json toJSON() const;
	bool operator==(const FieldReconstruction& rhs) const;
};

struct InCellPotentials {
	InCellPotentials() = default;
    InCellPotentials(const nlohmann::json&);

	bool operator==(const InCellPotentials&) const;

    Box innerRegionBox;
    std::map<MaterialId, FieldReconstruction> electric, magnetic;

    double getCapacitanceUsingInnerRegion(int i, int j) const;
    double getInductanceUsingInnerRegion(int i, int j) const;
    double getCapacitanceOnBox(int i, int j, const Box& box) const;
    double getInductanceOnBox(int i, int j, const Box& box) const;

    nlohmann::json toJSON() const;
};

}