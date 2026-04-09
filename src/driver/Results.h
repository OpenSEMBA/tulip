#pragma once

#include <nlohmann/json.hpp>

#include "FES.h"
#include "Model.h"
#include "multipolarExpansion.h"

namespace tulip {

void saveToJSONFile(const nlohmann::json&, const std::string& filename);

class MultiwireParameters {
public:
    virtual ~MultiwireParameters() = default;
    virtual nlohmann::json toJSON() const = 0;

    const Domain& getDomain() const { return domain; }
    void setDomain(const Domain& value) { domain = value; }
private:
    Domain domain;
};

class PULParameters : public MultiwireParameters {
public:
    PULParameters() = default;
    PULParameters(const nlohmann::json&);

    bool operator==(const PULParameters&) const;

    mfem::DenseMatrix getCapacitiveCouplingCoefficients() const;

    nlohmann::json toJSON() const;

    mfem::DenseMatrix L, C; // Stored in SI units.
    mfem::Vector R; // Stored in SI units.
};

class FieldReconstruction {
public:
    double innerRegionAveragePotential;
    std::array<double,2> expansionCenter;
    multipolarCoefficients ab; // Stored in natural units.
    std::map<ConductorId, double> conductorPotentials;

    nlohmann::json toJSON() const;
	bool operator==(const FieldReconstruction& rhs) const;
};

class InCellPotentials : public MultiwireParameters {
public:
	InCellPotentials() = default;
    InCellPotentials(const nlohmann::json&);

	bool operator==(const InCellPotentials&) const;

    double getInCellCapacitanceUsingInnerRegion(int i, int j) const;
    double getInCellInductanceUsingInnerRegion(int i, int j) const;
    double getInCellCapacitanceOnBox(int i, int j, const Box& box) const;
    double getInCellInductanceOnBox(int i, int j, const Box& box) const;

    const mfem::Array<double>& getResistancesPerMeter() const { return R; }
    mfem::Array<double>& getResistancesPerMeter() { return R; }
    const Box& getInnerRegionBox() const { return innerRegionBox; }
    Box& getInnerRegionBox() { return innerRegionBox; }
    const std::map<ConductorId, FieldReconstruction>& getElectric() const { return electric; }
    std::map<ConductorId, FieldReconstruction>& getElectric() { return electric; }
    const std::map<ConductorId, FieldReconstruction>& getMagnetic() const { return magnetic; }
    std::map<ConductorId, FieldReconstruction>& getMagnetic() { return magnetic; }

    nlohmann::json toJSON() const;

private:
    mfem::Array<double> R;
    Box innerRegionBox;
    std::map<ConductorId, FieldReconstruction> electric, magnetic;
};

class MultiwireParametersByDomain {
public:
    const InCellPotentials* getInCellPotentials() const;
    std::map<Domain::Id, const PULParameters*> getPULParameters() const;
    void setDomainTree(const DomainTree& value) { domainTree = value; }
    void add(Domain::Id id, std::unique_ptr<MultiwireParameters> value);
private:
    DomainTree domainTree;
    std::map<Domain::Id, std::unique_ptr<MultiwireParameters>> domainToPUL;
};

}