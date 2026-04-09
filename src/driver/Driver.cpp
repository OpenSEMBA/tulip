#include "Driver.h"

#include "Solver.h"
#include "AdaptedInputParser.h"
#include "multipolarExpansion.h"

using namespace mfem;

namespace tulip {

namespace {

bool ignoresDielectrics(Driver::FieldType fieldType)
{
	return fieldType == Driver::FieldType::magnetic;
}

const char* getFieldTypeName(Driver::FieldType fieldType)
{
	switch (fieldType) {
	case Driver::FieldType::electric:
		return "electric";
	case Driver::FieldType::magnetic:
		return "magnetic";
	}

	throw std::runtime_error("Unsupported field type.");
}

std::string getFieldTypeSuffix(Driver::FieldType fieldType)
{
	switch (fieldType) {
	case Driver::FieldType::electric:
		return "_electrostatic";
	case Driver::FieldType::magnetic:
		return "_magnetostatic";
	}

	throw std::runtime_error("Unsupported field type.");
}

}

void exportFieldSolutions(
	const DriverOptions& opts,
	Solver& s,
	int conductorId,
	Driver::FieldType fieldType,
	std::string extraName = "")
{
	const std::string suffix{ getFieldTypeSuffix(fieldType) };

	std::stringstream ss;
	ss << "Conductor_" << conductorId;

	if (opts.exportParaViewSolution) {
		std::string outputName{ opts.exportFolder + "/" + "ParaView/" + ss.str() + suffix + extraName};
		ParaViewDataCollection pd{ outputName, s.getMesh() };
		s.writeParaViewFields(pd);
	}
}


Driver Driver::loadFromAdaptedFile(const std::string& fn)
{
	AdaptedInputParser p{ fn };
	return Driver{
		p.readModel(),
		p.readDriverOptions()
	};
}

Driver::Driver(Model&& model, const DriverOptions& opts) :
	model_{ std::move(model) },
	opts_{ opts }
{
	if (model_.getMaterials().getConductors().empty()) { 
		throw std::runtime_error("Model must have at least one conductor.");
	}

	// Solve for all conductors.
	std::cout << "Solving electrostatic problems:" << std::endl;
	electric_ = solveForAllConductors(FieldType::electric);

	if (model_.getMaterials().hasDielectrics()) {
		std::cout << "Solving magnetostatic problems:" << std::endl;
		magnetic_ = solveForAllConductors(FieldType::magnetic);
	}
	else {
		std::cout << "No dielectrics found. Reusing electrostatic solution for magnetostatic." << std::endl;
	}
}

mfem::DenseMatrix Driver::getCFromGeneralizedC(
	const mfem::DenseMatrix& gC,
	const Model::Openness& opennness)
{
	// Implements (5.21) from,
	// Clayton Paul. Analysis of multiconductor transmision lines. 2007.
	mfem::DenseMatrix C(gC.NumRows() - 1, gC.NumCols() - 1);
	auto conductors = model_.getMaterials().getConductors();

	if (opennness == Model::Openness::closed) {
		for (int i = 1; i < gC.NumRows(); ++i) {
			for (int j = 1; j < gC.NumCols(); ++j) {
				C(i - 1, j - 1) = gC(i, j);
			}
		}
		return C;
	}

	double den = 0.0;
	for (int i = 0; i < gC.NumRows(); ++i) {
		for (int j = 0; j < gC.NumCols(); ++j) {
			den += gC(i, j);
		}
	}

	for (int i = 1; i <= C.NumRows(); ++i) {
		for (int j = 1; j <= C.NumCols(); ++j) {
			double rowSum = 0.0;
			for (int k = 0; k < gC.NumCols(); ++k) {
				rowSum += gC(i, k);
			}
			double colSum = 0.0;
			for (int m = 0; m < gC.NumRows(); ++m) {
				colSum += gC(m, j);
			}
			C(i - 1, j - 1) = gC(i, j) - rowSum * colSum / den;
		}
	}

	return C;
}

SolvedProblem Driver::solveForAllConductors(FieldType fieldType)
{
	SolvedProblem res;
	const auto baseParameters{ 
		buildSolverInputsFromModel(model_, fieldType) };
	res.solver = std::make_unique<Solver>(
		*model_.getMesh(), baseParameters, opts_.solverOptions);
	Solver& s = *res.solver.get();

	auto conductors = model_.getMaterials().getConductors();
	for (const auto c : conductors) {
		std::cout << "- Conductor #" 
			<< c->getConductorId() << "... " << std::flush;

		auto dbcs = baseParameters.dirichletBoundaries;
		dbcs[c->getAttribute()] = 1.0;
		s.setDirichletConditions(dbcs);
		s.Solve();

		exportFieldSolutions(opts_, s, c->getConductorId(), fieldType);
		res.solutions[c->getConductorId()] = std::move(s.getSolution());

		std::cout << "[OK]" << std::endl;
	}

	return res;
}

SolverInputs Driver::buildSolverInputsFromModel(
	const Model& model,
	FieldType fieldType)
{
	SolverInputs res;
	auto dielectrics = model.getMaterials().getDielectrics();
	for (const auto& d : dielectrics) {
		if (ignoresDielectrics(fieldType)) {
			res.domainPermittivities[d->getAttribute()] = 1.0;
		} else {
			res.domainPermittivities[d->getAttribute()] = d->getRelativePermittivity();
		}
	}

	for (auto b: model.getMaterials().getOpenBoundaries()) {
		res.openBoundaries.push_back(b->getAttribute());
	}
	
	for (auto b: model.getMaterials().getConductors()) {
		res.dirichletBoundaries[b->getAttribute()] = 0.0;
	}

	return res;
}

SolvedProblem* Driver::getSolvedProblem(FieldType fieldType)
{
	if (fieldType == FieldType::magnetic && model_.getMaterials().hasDielectrics()) {
		return &magnetic_;
	}

	return &electric_;
}

DenseMatrix Driver::getGeneralizedCMatrix(FieldType fieldType)
{
	// PUL generalized capacitance matrix for a N conductors system as defined in:
	// "Clayton Paul's book: Analysis of Multiconductor Transmission Lines"
	// - Generalized C contains N x N entries.

	// Preconditions. 
	const auto conductors{ model_.getMaterials().getConductors() };
	const auto openness{ model_.determineOpenness() };
	if (conductors.size() == 1 && openness == Model::Openness::closed) {
		throw std::runtime_error(
			"The number of conductors must be at least 2 for closed problems.");
	}

	int CSize = (int)conductors.size();
	mfem::DenseMatrix C(CSize);

	SolvedProblem* sP = getSolvedProblem(fieldType);

	int condI = 0;
	for (const auto& c : conductors) {
		sP->solver->setSolution(sP->solutions.at(c->getConductorId()));

		// Fills row
		int condJ = 0;
		for (const auto& c : conductors) {
			// C_ij = Q_j / V_i. V_i is always 1.0
			C(condI, condJ) = 
				sP->solver->getChargeInBoundary(c->getAttribute());
			condJ++;
		}

		exportFieldSolutions(
			opts_, *sP->solver, c->getConductorId(), fieldType);
		
			condI++;
	}

	C.Symmetrize();

	return C;
}

DenseMatrix Driver::getCMatrix()
{
	// PUL capacitance matrix for a N conductors system as defined in:
	// "Clayton Paul's book: Analysis of Multiconductor Transmission Lines"
	// - Standard C contains N-1 x N-1 entries
	
	auto gC = getGeneralizedCMatrix(FieldType::electric);
	auto res{ getCFromGeneralizedC(gC, model_.determineOpenness()) };
	return res;
}

DenseMatrix Driver::getLMatrix()
{
	// PUL inductance matrix as defined in:
	// Clayton Paul's book: Analysis of Multiconductor Transmission Lines
	// Contains N-1 x N-1 entries for a problem of N conductors.
	// Inductance matrix can be computed from the 
	// capacitance obtained ignoring dielectrics as
	//          L = mu0 * eps0 * C^{-1}
	auto gC = getGeneralizedCMatrix(FieldType::magnetic);
	auto res{ getCFromGeneralizedC(gC, model_.determineOpenness()) };
	res.Invert();
	return res;
}

PULParameters Driver::buildPULParametersForModel()
{
	PULParameters res;

	res.C = getCMatrix();
	res.C *= EPSILON0_SI;

	res.L = getLMatrix();
	res.L *= MU0_SI;

	return res;
}

PULParameters Driver::buildGeneralizedLCMatrices()
{
	PULParameters res;

	res.C = getGeneralizedCMatrix(FieldType::electric);
	res.C *= EPSILON0_SI;

	res.L = getGeneralizedCMatrix(FieldType::magnetic);
	res.L *= MU0_SI;

	return res;
}

void Driver::run()
{
	
	saveToJSONFile(
		getMultiwireParametersByDomains().toFDTDJSON(),
		opts_.exportFolder + "tulip.out.json");


}

PULParameters Driver::getPULMTL()
{	
	if (model_.getDomains().size() == 1) {
			return buildPULParametersForModel();
	} else {
		throw std::runtime_error("getPULMTL can only be called for single domain problems.");
	}
}

MultiwireParametersByDomain Driver::getMultiwireParametersByDomains()
{
	MultiwireParametersByDomain res;

	auto idToDomain{ model_.getDomains() };
	DomainTree domainTree{ idToDomain };
	res.setDomainTree(domainTree);
	const auto openness = model_.determineOpenness();
	std::unique_ptr<InCellPotentials> openDomainPotentials;
	if (openness == Model::Openness::open) {
		openDomainPotentials = std::make_unique<InCellPotentials>(getInCellPotentials());
	}

	// Build mapping from conductor ID to global matrix index.
	auto conductors = model_.getMaterials().getConductors();
	std::map<ConductorId, int> condIdToIndex;
	int idx = 0;
	for (const auto& c : conductors) {
		condIdToIndex[c->getConductorId()] = idx++;
	}

	// Get global generalized C matrices.
	auto globalGC = getGeneralizedCMatrix(FieldType::electric);
	auto globalGC0 = getGeneralizedCMatrix(FieldType::magnetic);

	for (const auto& [domId, domain] : idToDomain) {
		if (openness == Model::Openness::open &&
			domain.ground == Domain::UNDEFINED_GROUND) {
			auto inCell = std::make_unique<InCellPotentials>(*openDomainPotentials);

			auto restrictToDomain = [&](std::map<ConductorId, FieldReconstruction>& fields) {
				for (auto it = fields.begin(); it != fields.end();) {
					if (!domain.conductorIds.count(it->first)) {
						it = fields.erase(it);
						continue;
					}

					auto& potentials = it->second.conductorPotentials;
					for (auto pIt = potentials.begin(); pIt != potentials.end();) {
						if (!domain.conductorIds.count(pIt->first)) {
							pIt = potentials.erase(pIt);
						}
						else {
							++pIt;
						}
					}

					++it;
				}
			};

			restrictToDomain(inCell->electric);
			restrictToDomain(inCell->magnetic);
			inCell->setDomain(domain);
			res.add(domId, std::move(inCell));
			continue;
		}

		// Order conductors: ground first, then rest sorted.
		std::vector<ConductorId> orderedConds;
		ConductorId groundCond = (domain.ground != Domain::UNDEFINED_GROUND)
			? domain.ground
			: *domain.conductorIds.begin();
		orderedConds.push_back(groundCond);
		for (auto cId : domain.conductorIds) {
			if (cId != groundCond) {
				orderedConds.push_back(cId);
			}
		}

		int n = (int)orderedConds.size();

		// Build domain-local generalized C by aggregating interior conductors.
		// When exciting conductor i, set V=1 on i and all conductors inside it.
		// The effective charge on conductor j is the sum of charges on j and
		// all conductors inside it.
		auto buildDomainGC = [&](const mfem::DenseMatrix& gC) {
			mfem::DenseMatrix domGC(n);
			for (int i = 0; i < n; i++) {
				auto insideI = domainTree.getConductorsInsideConductor(orderedConds[i]);
				for (int j = 0; j < n; j++) {
					auto insideJ = domainTree.getConductorsInsideConductor(orderedConds[j]);
					double val = 0.0;
					for (auto m : insideI) {
						auto mIt = condIdToIndex.find(m);
						if (mIt == condIdToIndex.end()) continue;
						for (auto k : insideJ) {
							auto kIt = condIdToIndex.find(k);
							if (kIt == condIdToIndex.end()) continue;
							val += gC(mIt->second, kIt->second);
						}
					}
					domGC(i, j) = val;
				}
			}
			return domGC;
		};

		// Extract standard C (remove ground row/col).
		auto extractStdC = [&](const mfem::DenseMatrix& domGC) {
			mfem::DenseMatrix C(n - 1, n - 1);
			for (int i = 1; i < n; i++) {
				for (int j = 1; j < n; j++) {
					C(i - 1, j - 1) = domGC(i, j);
				}
			}
			return C;
		};

		auto pul = std::make_unique<PULParameters>();
		pul->setDomain(domain);

		pul->C = extractStdC(buildDomainGC(globalGC));
		pul->C *= EPSILON0_SI;

		pul->L = extractStdC(buildDomainGC(globalGC0));
		pul->L.Invert();
		pul->L *= MU0_SI;

		res.add(domId, std::move(pul));
	}

	return res;
}

std::map<ConductorId, double> Driver::getFloatingPotentials(
	ConductorId prescribedId, FieldType fieldType)
{
	// Returns a map from conductorId to value, with the potentials
	// which other conductors must have in order to be "floating",
	// i.e. with zero charge.
	// For open problems returns a map with N entries.

	auto conductors = model_.getMaterials().getConductors();

	bool prescribedIdExists = false;
	for (auto c : conductors ) {
		if (c->getConductorId() == prescribedId) {
			prescribedIdExists = true;
			break;
		}
	}
	if (!prescribedIdExists) {
		throw std::runtime_error("Invalid prescribed id");
	}
	
	if (model_.determineOpenness() == Model::Openness::closed) {
		throw std::runtime_error("Floating potentials not implemented for open problems");
	}


	// Special case for one conductor.
	if (conductors.size() == 1) {
		std::map<ConductorId, double> res;
		auto cId = conductors.front()->getConductorId();
		res[cId] = 1.0;
		return res;
	}

	// Determine C matrix.
	mfem::DenseMatrix C;
	auto openness = model_.determineOpenness();
	switch (openness) {
	case Model::Openness::closed:
		C = getCFromGeneralizedC(
				getGeneralizedCMatrix(fieldType), openness);
		break;
	case Model::Openness::open:
		C = getGeneralizedCMatrix(fieldType);
		break;
	default:
		throw std::runtime_error(
			"Floating potentials not implemented for this kind of openness.");
	}

	return computeFloatingPotentialsFromC(prescribedId, C);
}

std::map<ConductorId, double> Driver::computeFloatingPotentialsFromC(
	ConductorId prescribedId,
	const mfem::DenseMatrix& C) const
{
	std::map<ConductorId, double> res;
	auto conductors = model_.getMaterials().getConductors();
	
	// Forms system of equations to determine floating potentials. 
	// Q1 = C11*V1 + C21*V2 + ....
	// For prescribed V_1 = 1.0 we have C V = Q
	//    [ C ] [1.0, V_2, ...]^T = [Q_1, 0.0, ...]
	// 
	// which can be converted to A x = b with unknowns x = [Q_1, V_2, ...]^T
	auto N = C.NumRows();
	
	mfem::DenseMatrix A{ C };
	mfem::Vector b(N);
	mfem::Vector x(N);
	mfem::Vector negativeQ(N);
	
	int iPrescribed = 0;
	int i = 0;
	for (auto cI : conductors) {
		if (cI->getConductorId() == prescribedId) {
			iPrescribed = i;
			break;
		}
		i++;
	}
	
	negativeQ = 0.0;
	negativeQ(iPrescribed) = -1.0;
	
	A.SetCol(iPrescribed, negativeQ);

	b = C.GetColumn(iPrescribed);
	b *= -1.0;

	mfem::DenseMatrixInverse Ainv(A);
	Ainv.Mult(b, x);

	i = 0;
	for (auto c: conductors) {
		if (i == iPrescribed) {
			res[c->getConductorId()] = 1.0;
		} else {
			res[c->getConductorId()] = x(i);
		}
		i++;
	}

	return res;
}

double Driver::getInnerRegionAveragePotential(
	const Solver& s,
	bool includeConductors)
{

	double totalPotential = 0.0;
	double totalArea = 0.0;
	
	for (const auto& m: model_.getMaterials().getAll() ) {
		if (model_.isOuterRegion(m)) {
			continue;
		}
		auto area = model_.getAreaOfMaterial(m); 
		if (m->isDomainMaterial()) {
			totalPotential += 
				s.getAveragePotentialInDomain(m->getAttribute()) * area;
		} else {
			totalPotential += 
				s.getAveragePotentialInBoundary(m->getAttribute()) * area;
		}
		totalArea += area;
	}

	return totalPotential / totalArea;
}

void Driver::loadFloatingPotentials(
	SolvedProblem* sP, const std::map<ConductorId, double>& fp) const
{
	Solver& s = *sP->solver;
	
	s.getPhi() *= 0.0;
	s.getE() *= 0.0;
	s.getD() *= 0.0;
	for (const auto& cJ : model_.getMaterials().getConductors()) {
		auto condJ = cJ->getConductorId();
		auto fpJ = fp.at(condJ);
		s.getPhi().Add(fpJ, *sP->solutions[condJ].phi);
		s.getE().Add(fpJ, *sP->solutions[condJ].e);
		s.getD().Add(fpJ, *sP->solutions[condJ].d);
	}
}

std::map<ConductorId, FieldReconstruction> Driver::getFieldParameters(
	FieldType fieldType)
{
	std::map<ConductorId, FieldReconstruction> res;

	SolvedProblem* sP = getSolvedProblem(fieldType);

	Solver& s = *sP->solver;

	std::cout << "- Computing " << getFieldTypeName(fieldType)
		<< " field coefficients." << std::endl;
	const auto conductors = model_.getMaterials().getConductors();

	// Compute the C matrix once for all conductors to avoid
	// reassembling operators N times (was O(N^2), now O(N)).
	mfem::DenseMatrix C = getGeneralizedCMatrix(fieldType);

	for (const auto& cI : conductors) {
		
		auto condI = cI->getConductorId();

		std::cout << "- Conductor #" << condI << "... " << std::flush;
		auto fp = computeFloatingPotentialsFromC(condI, C);
		
		loadFloatingPotentials(sP, fp);

		exportFieldSolutions(opts_, s, condI, fieldType, "_floating");

		res[condI].innerRegionAveragePotential = 
			getInnerRegionAveragePotential(s, true);
		auto centerOfCharge = s.getCenterOfCharge();
		std::copy(
			centerOfCharge.begin(), centerOfCharge.end(), 
			res[condI].expansionCenter.begin());
		res[condI].ab = s.getMultipolarCoefficients(opts_.multipolarExpansionOrder);
		for (const auto& cJ : conductors) {
			auto condJ = cJ->getConductorId();
			res[condI].conductorPotentials[condJ] = fp.at(condJ);
		}
		std::cout << "[OK]" << std::endl;
	}

	return res;
}

InCellPotentials Driver::getInCellPotentials()
{
	if (model_.determineOpenness() != Model::Openness::open) {
		throw std::runtime_error("In cell potentials can only be determined for open problems.");
	}

	InCellPotentials res;
	res.innerRegionBox = model_.getInnerRegionBoundingBox();

	res.electric = getFieldParameters(FieldType::electric);
	res.magnetic = getFieldParameters(FieldType::magnetic);

	return res;
}

}