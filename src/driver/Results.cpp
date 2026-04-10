#include "Results.h"

#include "constants.h"

#include <filesystem>
#include <fstream>

namespace tulip {

using namespace mfem;
using json = nlohmann::json;

std::vector<double> toVec(const Vector& vec)
{
	std::vector<double> r(vec.Size());
	for (auto i{ 0 }; i < vec.Size(); ++i) {
        r[i] = vec[i];
	}
	return r;
}

std::vector<std::vector<double>> toVecVec(const DenseMatrix& m)
{
    std::vector<std::vector<double>> r(m.NumRows());
    for (auto i{ 0 }; i < m.NumRows(); ++i) {
        r[i].resize(m.NumCols());
        for (auto j{ 0 }; j < m.NumCols(); ++j) {
            r[i][j] = m(i, j);
        }
    }

    return r;
}

Vector toMFEMVector(const std::vector<double>& v)
{
    Vector r(v.size());
	for (auto i{ 0 }; i < v.size(); ++i) {
        r[i] = v[i];
	}
	return r;
}

DenseMatrix toMFEMDenseMatrix(const std::vector<std::vector<double>>& v)
{
    
    if (v.size() == 0 || v[0].size() == 0) {
        return {};
    }

    for (const auto& row: v) {
        if (row.size() != v[0].size()) {
            throw std::runtime_error("vec-vec matrix is not square.");
        }
    }
    
    DenseMatrix r(int(v.size()), int(v[0].size()));
    for (auto i{ 0 }; i < r.NumRows(); ++i) {
        for (auto j{ 0 }; j < r.NumCols(); ++j) {
            r(i, j) = v[i][j];
        }
    }

    return r;
}

Vector zeroVectorLike(const DenseMatrix& matrix)
{
    Vector result(matrix.NumRows());
    for (auto i{ 0 }; i < result.Size(); ++i) {
        result[i] = 0.0;
    }
    return result;
}

PULParameters::PULParameters(const json& j)
{
    if (!j.contains("C") || !j.contains("L")) {
        throw std::runtime_error("JSON does not contain C and/or L.");
    }
    C = toMFEMDenseMatrix(j["C"]);
    L = toMFEMDenseMatrix(j["L"]);
    if (j.contains("R")) {
        if (j["R"].is_array() && !j["R"].empty() && j["R"][0].is_array()) {
            auto backwardCompatibleR = toMFEMDenseMatrix(j["R"]);
            R.SetSize(backwardCompatibleR.NumRows());
            for (int i = 0; i < backwardCompatibleR.NumRows(); ++i) {
                R[i] = backwardCompatibleR(i, i);
            }
        }
        else {
            R = toMFEMVector(j["R"]);
        }
    }
    else {
        R = zeroVectorLike(C);
    }
}

bool PULParameters::operator==(const PULParameters& rhs) const
{ 
    return
		toVec(R) == toVec(rhs.R) &&
        toVecVec(C) == toVecVec(rhs.C) &&
        toVecVec(L) == toVecVec(rhs.L);
}

DenseMatrix PULParameters::getCapacitiveCouplingCoefficients() const
{
    DenseMatrix r(C.NumRows(), C.NumCols());
    for (auto i{ 0 }; i < C.NumRows(); ++i) {
        auto selfC{ C(i,i) };
        for (auto j{ 0 }; j < C.NumCols(); ++j) {
            r(i, j) = C(i, j) / selfC;
        }
    }
    return r;
}

json PULParameters::toJSON() const
{
    json res;
    res["R"] = toVec(R);
    res["C"] = toVecVec(C);
    res["L"] = toVecVec(L);
    return res;
}

void saveToJSONFile(const json& j, const std::string& filename)
{
    const auto outputPath = std::filesystem::path(filename);
    if (outputPath.has_parent_path()) {
        std::filesystem::create_directories(outputPath.parent_path());
    }

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        throw std::runtime_error("Unable to open output file: " + filename);
    }
    ofs << j;
}

multipolarCoefficients multipolarCoefficientsFromJSON(const json& j)
{
    multipolarCoefficients res;
	for (const auto& coeffs : j) {
        res.push_back(std::make_pair(coeffs[0], coeffs[1]));
	}
	return res;
}

std::map<ConductorId, FieldReconstruction> potentialsFromJSON(const json& j)
{
    std::map<ConductorId, FieldReconstruction> res;
    int reconstructedId = 0;
    for (const auto& f : j) {
		FieldReconstruction fr;
		
        fr.innerRegionAveragePotential = f.at("innerRegionAveragePotential").get<double>();
		fr.expansionCenter = f.at("expansionCenter").get<std::array<double,2>>();
		
        fr.ab = multipolarCoefficientsFromJSON(f.at("ab"));
        
        int id = 0;
        for (const auto& potential : f.at("conductorPotentials")) {
            fr.conductorPotentials[id++] = potential;
		}

		res[reconstructedId++] = fr;
    }
    
    return res;
}

InCellPotentials::InCellPotentials(const nlohmann::json& j)
{
    innerRegionBox.min = j.at("innerRegionBox").at("min").get<std::array<double,2>>();
	innerRegionBox.max = j.at("innerRegionBox").at("max").get<std::array<double,2>>();

    electric = potentialsFromJSON(j.at("electric"));
	magnetic = potentialsFromJSON(j.at("magnetic"));
}

bool InCellPotentials::operator==(const InCellPotentials& rhs) const
{
    
    if (!(innerRegionBox == rhs.innerRegionBox)) {
        return false;
    }

    if (!(electric == rhs.electric)) {
        return false;
    }

    if (!(magnetic == rhs.magnetic)) {
        return false;
    }
    return true;
}

double InCellPotentials::getInCellCapacitanceUsingInnerRegion(int i, int j) const
{
    double Qj = electric.at(j).ab[0].first;
    double avVj = electric.at(j).innerRegionAveragePotential;
    double ViWhenPrescribedVj = electric.at(j).conductorPotentials.at(i);
    avVj = -avVj + ViWhenPrescribedVj;

    return Qj / avVj * EPSILON0_SI;
}

double InCellPotentials::getInCellInductanceUsingInnerRegion(int i, int j) const
{
    double Ij = magnetic.at(j).ab[0].first;
    double avAj = magnetic.at(j).innerRegionAveragePotential;
    double AiWhenPrescribedAj = magnetic.at(j).conductorPotentials.at(i);
    avAj = -avAj + AiWhenPrescribedAj;
    
    return avAj / Ij * MU0_SI;
}

std::array<std::vector<double>, 2> buildIntegrationPlanesForBox(
    const Box& integrationBox, 
    const Box& innerRegionBox) 
{
    const int GRID_INTEGRATION_SAMPLING_POINTS = 100;

    for (int x = 0; x < 2; ++x) {
        if ((integrationBox.min[x] > innerRegionBox.min[x]) ||
            (integrationBox.max[x] < innerRegionBox.max[x])) {
            throw std::runtime_error("Integration box has to be larger than the  inner region.");
        }
    }

    std::array<std::set<double>, 2> planes;
    for (int x = 0; x < 2; ++x) {
        std::set<double> controlPoints;
        controlPoints.insert(integrationBox.min[x]);
        controlPoints.insert(integrationBox.max[x]);
        controlPoints.insert(innerRegionBox.min[x]);
        controlPoints.insert(innerRegionBox.max[x]);

        // Fills with equispaced points
        for (auto it = std::next(controlPoints.begin()); it != controlPoints.end(); ++it) {
            auto prev = std::prev(it);
            const double step = (*it - *prev) / GRID_INTEGRATION_SAMPLING_POINTS;
            for (int k = 0; k < GRID_INTEGRATION_SAMPLING_POINTS; ++k) {
                const double newPoint = *prev + k * step;
                planes[x].insert(newPoint);
            }
			planes[x].insert(*it); 
        }
    }

	std::array<std::vector<double>, 2> res;
	for (int x = 0; x < 2; ++x) {
		res[x].assign(planes[x].begin(), planes[x].end());
	}
	return res;
}

double getAveragePotential(
    const FieldReconstruction& potential, 
    const Box& innerBox, 
    const Box& outerBox)
{
    auto integrationPlanes{ buildIntegrationPlanesForBox(outerBox, innerBox) };
    double outerV = 0.0;
    for (int m = 1; m < integrationPlanes[0].size(); ++m) {
        for (int n = 1; n < integrationPlanes[1].size(); ++n) {
            const auto& xMin = integrationPlanes[0][m - 1];
            const auto& xMax = integrationPlanes[0][m];
            const auto& yMin = integrationPlanes[1][n - 1];
            const auto& yMax = integrationPlanes[1][n];

            mfem::Vector midPoint({ 0.5 * (xMin + xMax), 0.5 * (yMin + yMax) });
            double area = (xMax - xMin) * (yMax - yMin);
            if (innerBox.isWithinBox(midPoint)) {
                // If the point is within the inner region, we do not integrate over it.
                continue;
            }
            mfem::Vector center(2);
            std::copy(
                potential.expansionCenter.begin(), potential.expansionCenter.end(),
                center.begin());
            outerV += area * multipolarExpansion(midPoint, potential.ab, center);
        }
    }

    double innerV = potential.innerRegionAveragePotential * innerBox.area();
    double avVj = (innerV + outerV) / outerBox.area();

    return avVj;
}

double InCellPotentials::getInCellCapacitanceOnBox(int i, int j, const Box& cellBox) const
{
	double Qj = electric.at(j).ab[0].first;

    double avVj = getAveragePotential(electric.at(j), innerRegionBox, cellBox);
	double ViWhenPrescribedVj = electric.at(j).conductorPotentials.at(i);
	avVj = -avVj + ViWhenPrescribedVj;
	
	return Qj / avVj * EPSILON0_SI ;
}

double InCellPotentials::getInCellInductanceOnBox(int i, int j, const Box& cellBox) const
{
    double Ij = magnetic.at(j).ab[0].first;

    double avAj = getAveragePotential(magnetic.at(j), innerRegionBox, cellBox);
    double AiWhenPrescribedAj = magnetic.at(j).conductorPotentials.at(i);
    avAj = -avAj + AiWhenPrescribedAj;

    return avAj / Ij * MU0_SI;
}

nlohmann::json InCellPotentials::toJSON() const
{
    nlohmann::json res;
    res["innerRegionBox"] = nlohmann::json::object({
        {"min", innerRegionBox.min},
        {"max", innerRegionBox.max}
    });

    res["electric"] = nlohmann::json::array();
    for (const auto& [matId, fieldReconstruction] : electric) {
		
        res["electric"].push_back(electric.at(matId).toJSON());
	}
	
    res["magnetic"] = nlohmann::json::array();
	for (const auto& [matId, fieldReconstruction] : magnetic) {
		res["magnetic"].push_back(magnetic.at(matId).toJSON());
	}

    return res;
}

nlohmann::json FieldReconstruction::toJSON() const
{
	nlohmann::json res;

	res["innerRegionAveragePotential"] = innerRegionAveragePotential;
	res["expansionCenter"] = expansionCenter;
	res["ab"] = nlohmann::json::array();
	for (const auto& [a, b] : ab) {
		res["ab"].push_back({ a, b });
	}
    res["conductorPotentials"] = nlohmann::json::array();
    for (const auto& [condId, v] : conductorPotentials) {
        res["conductorPotentials"].push_back(v);
    }

    return res;
}

bool FieldReconstruction::operator==(const FieldReconstruction& rhs) const
{
    bool res;

	res = innerRegionAveragePotential == rhs.innerRegionAveragePotential;
    
    res &= expansionCenter[0] == rhs.expansionCenter[0];
	res &= expansionCenter[1] == rhs.expansionCenter[1];

	res &= ab == rhs.ab;
	
    if (conductorPotentials.size() != rhs.conductorPotentials.size()) {
		return false;
	}
	for (const auto& [condId, v] : conductorPotentials) {
		if (rhs.conductorPotentials.find(condId) == rhs.conductorPotentials.end()) {
			return false;
		}
		res &= v == rhs.conductorPotentials.at(condId);
	}
	return res;
}

const InCellPotentials* MultiwireParametersByDomain::getInCellPotentials() const
{
    for (const auto& [domainId, params] : domainToPUL) {
        (void)domainId;
        if (auto* inCell = dynamic_cast<InCellPotentials*>(params.get())) {
            return inCell;
        }
    }
    return nullptr;
}

std::map<Domain::Id, const PULParameters*> MultiwireParametersByDomain::getPULParameters() const
{
    std::map<Domain::Id, const PULParameters*> res;
    for (const auto& [domainId, params] : domainToPUL) {
        if (auto* pul = dynamic_cast<PULParameters*>(params.get())) {
            res[domainId] = pul;
        }
    }
    return res;
}

void MultiwireParametersByDomain::add(
    Domain::Id id,
    std::unique_ptr<MultiwireParameters> value)
{
    domainToPUL[id] = std::move(value);
}

nlohmann::json MultiwireParametersByDomain::toFDTDJSON() const
{
    json result;
    json materials = json::array();
    json materialAssociations = json::array();

    int materialId = 1;
    for (const auto& [domainId, params] : domainToPUL) {
        json mat;
        mat["id"] = materialId;

        if (auto* pul = dynamic_cast<PULParameters*>(params.get())) {
            auto pulJSON = pul->toJSON();
            mat["type"] = "shieldedMultiwire";
            mat["capacitancePerMeter"] = pulJSON["C"];
            mat["inductancePerMeter"] = pulJSON["L"];
            mat["resistancePerMeter"] = pulJSON["R"];
        }
        else if (auto* inCell = dynamic_cast<InCellPotentials*>(params.get())) {
            mat["type"] = "unshieldedMultiwire";
            mat["inCellParameters"]["multipolarExpansion"] = inCell->toJSON();
        }

        materials.push_back(mat);

        json assoc;
        assoc["materialId"] = materialId;
        json elementIds = json::array();
        const auto& domain = params->getDomain();
        if (dynamic_cast<PULParameters*>(params.get()) != nullptr) {
            const auto groundCond = (domain.ground != Domain::UNDEFINED_GROUND)
                ? domain.ground
                : *domain.conductorIds.begin();
            for (auto condId : domain.conductorIds) {
                if (condId != groundCond) {
                    elementIds.push_back(condId);
                }
            }
        }
        else {
            for (auto condId : domain.conductorIds) {
                elementIds.push_back(condId);
            }
        }
        assoc["elementIds"] = elementIds;
        if (dynamic_cast<PULParameters*>(params.get()) != nullptr &&
            domain.ground != Domain::UNDEFINED_GROUND) {
            assoc["containedWithinElementId"] = domain.ground;
        }
        materialAssociations.push_back(assoc);

        materialId++;
    }

    result["materials"] = materials;
    result["materialAssociations"] = materialAssociations;
    return result;
}

}