#include "Driver.h"

#include "Solver.h"
#include "AdaptedInputParser.h"
#include "multipolarExpansion.h"

using namespace mfem;

namespace tulip {

void exportFieldSolutions(
	const DriverOptions& opts,
	Solver& s,
	int conductorId,
	bool ignoreDielectrics,
	std::string extraName = "")
{
	const std::string suffix{ ignoreDielectrics ? "_magnetostatic" : "_electrostatic" };

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
	electric_ = solveForAllConductors(false);

	if (model_.getMaterials().hasDielectrics()) {
		std::cout << "Solving magnetostatic problems." << std::endl;
		magnetic_ = solveForAllConductors(true);
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
		int i = 0;
		auto groundId = conductors.front()->getConductorId();
		for (auto cI : conductors) {
			if (cI->getConductorId() == groundId) {
				continue;
			}
			int j = 0;
			for (auto cJ : conductors) {
				if (cJ->getConductorId() == groundId) {
					continue;
				}
				C(i - 1, j - 1) = gC(i, j);
				j++;
			}
			i++;
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

SolvedProblem Driver::solveForAllConductors(bool ignoreDielectrics)
{
	SolvedProblem res;
	const auto baseParameters{ 
		buildSolverInputsFromModel(model_, ignoreDielectrics) };
	res.solver = std::make_unique<Solver>(
		*model_.getMesh(), baseParameters, opts_.solverOptions);
	Solver& s = *res.solver.get();

	auto conductors = model_.getMaterials().getConductors();
	for (const auto c : conductors) {
		std::cout << "- Solving conductor #" 
			<< c->getConductorId() << "... " << std::flush;

		auto dbcs = baseParameters.dirichletBoundaries;
		dbcs[c->getAttribute()] = 1.0;
		s.setDirichletConditions(dbcs);
		s.Solve();

		exportFieldSolutions(opts_, s, c->getConductorId(), ignoreDielectrics);
		res.solutions[c->getConductorId()] = std::move(s.getSolution());

		std::cout << "[OK]" << std::endl;
	}

	return res;
}

SolverInputs Driver::buildSolverInputsFromModel(
	const Model& model,
	bool ignoreDielectrics)
{
	SolverInputs res;
	auto dielectrics = model.getMaterials().getDielectrics();
	for (const auto& d : dielectrics) {
		if (ignoreDielectrics) {
			res.domainPermittivities.at(d->getAttribute()) = 1.0;
		} else {
			res.domainPermittivities.at(d->getAttribute()) = d->getRelativePermittivity();
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

DenseMatrix Driver::getGeneralizedCMatrix(bool ignoreDielectrics)
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

	SolvedProblem* sP;
	if (ignoreDielectrics) {
		sP = model_.getMaterials().hasDielectrics() ? &magnetic_ : &electric_;
	}
	else {
		sP = &electric_;
	}

	int condI = 0;
	for (const auto& c : conductors) {
		sP->solver->setSolution(sP->solutions[condI]);

		// Fills row
		int condJ = 0;
		for (const auto& c : conductors) {
			// C_ij = Q_j / V_i. V_i is always 1.0
			C(condI, condJ) = 
				sP->solver->getChargeInBoundary(c->getAttribute());
			condJ++;
		}

		exportFieldSolutions(
			opts_, *sP->solver, c->getConductorId(), ignoreDielectrics);
		
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
	
	auto gC = getGeneralizedCMatrix(false);
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
	auto gC = getGeneralizedCMatrix(true);
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

	res.C = getGeneralizedCMatrix(false);
	res.C *= EPSILON0_SI;

	res.L = getGeneralizedCMatrix(true);
	res.L *= MU0_SI;

	return res;
}

void Driver::run()
{
	auto openness{ model_.determineOpenness() };
	if (openness == Model::Openness::closed) {
		auto pul = buildPULParametersForModel();
		saveToJSONFile(
			pul.toJSON(), 
			opts_.exportFolder + "tulip.out.json");
	}
	else if (openness == Model::Openness::open) {
		auto inCell = getInCellPotentials();
		saveToJSONFile(
			inCell.toJSON(),
			opts_.exportFolder + "inCellPotentials.out.json");

		auto generalizedLCMatrices = buildGeneralizedLCMatrices();
		saveToJSONFile(
			generalizedLCMatrices.toJSON(),
			opts_.exportFolder + "generalizedLC.out.json");
	}
	else {
		throw std::runtime_error("Openness of the model is not supported.");
	}
}

PULParameters Driver::getPULMTL()
{
	return buildPULParametersForModel();
}

PULParametersByDomain Driver::getPULMTLByDomains()
{
	PULParametersByDomain res;

	auto idToDomain{ Domain::buildDomains(model_) };

	for (const auto& [id, domain] : idToDomain) {
		auto globalMesh{ *model_.getMesh() };
		auto domainModel = Domain::buildModelForDomain(globalMesh, model_.getMaterials(), domain);
		Driver subDomainDriver(std::move(domainModel),opts_);
		res.domainToPUL[id] = subDomainDriver.getPULMTL();
	}

	res.domainTree = DomainTree{ idToDomain };

	return res;
}

std::map<ConductorId, double> Driver::getFloatingPotentials(
	ConductorId prescribedId, bool ignoreDielectrics)
{
	// Returns a map from conductorId to value, with the potentials
	// which other conductors must have in order to be "floating",
	// i.e. with zero charge.
	// For open problems returns a map with N entries.

	std::map<ConductorId, double> res;

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
				getGeneralizedCMatrix(ignoreDielectrics), openness);
		break;
	case Model::Openness::open:
		C = getGeneralizedCMatrix(ignoreDielectrics);
		break;
	default:
		throw std::runtime_error(
			"Floating potentials not implemented for this kind of openness.");
	}
	
	// Forms system of equations to determine floating potentials. 
	// Q1 = C11*V1 + C21*V2 + ....
	// For prescribed V_1 = 1.0 we have C V = Q
	//    [ C ] [1.0, V_2, ...]^T = [Q_1, 0.0, ...]
	// 
	// which can be converted to A x = b with unknowns x = [Q_1, V_2, ...]^T
	auto N = C.NumRows();
	
	mfem::DenseMatrix A{ C };
	mfem::Vector negativeQ(N);
	negativeQ = 0.0;

	int i = 0;
	for (auto cI : conductors) {
		if (cI->getConductorId() == prescribedId) {
			negativeQ(i) = -1.0;
		}
		i++;
	}

	A.SetCol(i, negativeQ);

	mfem::Vector b(N);
	b = C.GetColumn(i);
	b *= -1.0;

	mfem::Vector x(N);
	mfem::DenseMatrixInverse Ainv(A);
	Ainv.Mult(b, x);

	i = 0;
	for (auto c: conductors) {
		res[c->getConductorId()] = x(i);
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
		if (m->isOuterRegion()) {
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

std::map<ConductorId, FieldReconstruction> Driver::getFieldParameters(
	bool ignoreDielectrics)
{
	std::map<ConductorId, FieldReconstruction> res;

	SolvedProblem* sP;
	if (ignoreDielectrics) {
		sP = model_.getMaterials().hasDielectrics() ? &magnetic_ : &electric_;
	}
	else {
		sP = &electric_;
	}

	Solver& s = *sP->solver;

	const auto conductors = model_.getMaterials().getConductors();
	for (const auto& cI : conductors) {
		
		auto condI = cI->getConductorId();
		auto fp = getFloatingPotentials(condI, ignoreDielectrics);
		s.getPhi() *= 0.0;
		s.getE() *= 0.0;
		s.getD() *= 0.0;
		for (const auto& cJ : conductors) {
			auto condJ = cJ->getConductorId();
			auto fpJ = fp.at(condJ);
			s.getPhi().Add(fpJ, *sP->solutions[condJ].phi);
			s.getE().Add(fpJ, *sP->solutions[condJ].e);
			s.getD().Add(fpJ, *sP->solutions[condJ].d);
		}

		exportFieldSolutions(opts_, s, condI, ignoreDielectrics, "_floating");

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
	}

	return res;
}

InCellPotentials Driver::getInCellPotentials()
{
	InCellPotentials res;

	if (model_.determineOpenness() != Model::Openness::open) {
		throw std::runtime_error("In cell parameters can only be computed for open problems.");
	}

	res.innerRegionBox = model_.getInnerRegionBoundingBox();

	res.electric = getFieldParameters(false);
	res.magnetic = getFieldParameters(true);

	return res;
}

}