#include "Driver.h"

#include "ElectrostaticSolver.h"
#include "Parser.h"
#include "multipolarExpansion.h"

using namespace mfem;

namespace pulmtln {

const std::string INNER_REGION_DEFAULT_NAME = "Vacuum_0";

AttrToValueMap buildAttrToValueMap(
	const NameToAttrMap& matToAtt, double value)
{
	AttrToValueMap bcs;
	for (const auto& [mat, att] : matToAtt) {
		bcs[att] = value;
	}
	return bcs;
}

std::vector<int> getAttributesInMap(const NameToAttrMap& m)
{
	std::vector<int> res;
	for (const auto& [k, v] : m) {
		res.push_back(v);
	}
	return res;
}

void exportFieldSolutions(
	const DriverOptions& opts,
	ElectrostaticSolver& s,
	const std::string name,
	bool ignoreDielectrics)
{
	if (opts.exportParaViewSolution) {
		std::string outputName{ opts.exportFolder + "/" + "ParaView/" + name };
		if (ignoreDielectrics) {
			outputName += "_no_dielectrics";
		}
		ParaViewDataCollection pd{ outputName, s.getMesh() };
		s.writeParaViewFields(pd);
	}

	if (opts.exportVisItSolution) {
		std::string outputName{ opts.exportFolder + "/" + "VisIt/" + name };
		if (ignoreDielectrics) {
			outputName += "_no_dielectrics";
		}
		VisItDataCollection dC{ outputName, s.getMesh() };
		s.writeVisItFields(dC);
	}
}


Driver Driver::loadFromFile(const std::string& fn)
{
	Parser p{ fn };
	return Driver{
		p.readModel(),
		p.readDriverOptions()
	};
}

Driver::Driver(Model&& model, const DriverOptions& opts) :
	model_{ std::move(model) },
	opts_{ opts }
{
	auto conductors{ model_.getMaterials().buildNameToAttrMapFor<PEC>() };
	std::vector<int> conductorIds;
	conductorIds.reserve(conductors.size());
	for (const auto& [name, attr] : conductors) {
		conductorIds.push_back(Materials::getMaterialIdFromName(name));
	}
	std::sort(conductorIds.begin(), conductorIds.end());

	// Preconditions.
	if (conductorIds.empty()) { 
		throw std::runtime_error("Model must have at least one conductor.");
	}

	if (conductorIds.front() != 0) {
		throw std::runtime_error(
			"Conductor with id 0 must be present in the model. ");
	}

	for (int i = 1; i < conductorIds.size(); ++i) {
		if (conductorIds[i] != conductorIds[i - 1] + 1) {
			throw std::runtime_error(
				"Conductor ids must be consecutive.");
		}
	}

	// Solve for all conductors.
	std::cout << "Solving electrostatic problems:" << std::endl;
	electric_ = solveForAllConductors(false);

	std::cout << "Solving magnetostatic problems." << std::endl;
	magnetic_ = solveForAllConductors(true);
}

mfem::DenseMatrix Driver::getCFromGeneralizedC(
	const mfem::DenseMatrix& gC,
	const Model::Openness& opennness)
{
	// Implements (5.21) from,
	// Clayton Paul. Analysis of multiconductor transmision lines. 2007.
	
	
	mfem::DenseMatrix C(gC.NumRows() - 1, gC.NumCols() - 1);

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

SolvedProblem Driver::solveForAllConductors(bool ignoreDielectrics)
{
	SolvedProblem res;
	const auto baseParameters{ 
		buildSolverInputsFromModel(model_, ignoreDielectrics) };
	res.solver = std::make_unique<ElectrostaticSolver>(
		*model_.getMesh(), baseParameters, opts_.solverOptions);
	ElectrostaticSolver& s = *res.solver.get();

	auto conductors{ model_.getMaterials().buildIdToAttrMapFor<PEC>() };
	res.solutions.resize(conductors.size());
	for (const auto& [condI, bdrAttI] : conductors) {
		std::cout << "- Solving conductor #" << condI << "... " << std::flush;

		auto dbcs = baseParameters.dirichletBoundaries;
		dbcs[bdrAttI] = 1.0;
		s.setDirichletConditions(dbcs);
		s.Solve();

		exportFieldSolutions(opts_, s, model_.getMaterials().get<PEC>(condI).name, ignoreDielectrics);
		res.solutions[condI] = std::move(s.getSolution());

		std::cout << "[OK]" << std::endl;
	}

	return res;
}

SolverInputs Driver::buildSolverInputsFromModel(
	const Model& model,
	bool ignoreDielectrics)
{
	const Materials& mats = model.getMaterials();

	SolverInputs parameters;

	auto domainToEpsr{
		buildAttrToValueMap(mats.buildNameToAttrMapFor<Dielectric>(), 1.0)
	};
	if (!ignoreDielectrics) {
		for (const auto& d : mats.dielectrics) {
			domainToEpsr.at(d.attribute) = d.relativePermittivity;
		}
	}
	parameters.domainPermittivities = domainToEpsr;
	parameters.openBoundaries = getAttributesInMap(mats.buildNameToAttrMapFor<OpenBoundary>());
	parameters.dirichletBoundaries = buildAttrToValueMap(mats.buildNameToAttrMapFor<PEC>(), 0.0);

	return parameters;
}

DenseMatrix Driver::getGeneralizedCMatrix(bool ignoreDielectrics)
{
	// PUL generalized capacitance matrix for a N conductors system as defined in:
	// "Clayton Paul's book: Analysis of Multiconductor Transmission Lines"
	// - Generalized C contains N x N entries.

	// Preconditions. 
	const auto conductors{ model_.getMaterials().buildIdToAttrMapFor<PEC>() };
	const auto openness{ model_.determineOpenness() };
	if (conductors.size() == 1 && openness == Model::Openness::closed) {
		throw std::runtime_error(
			"The number of conductors must be at least 2 for closed problems.");
	}

	int CSize = (int)conductors.size();
	mfem::DenseMatrix C(CSize);

	SolvedProblem* sP;
	if (ignoreDielectrics) {
		sP = &magnetic_;
	}
	else {
		sP = &electric_;
	}

	for (const auto& [condI, bdrAttI] : conductors) {
		sP->solver->setSolution(sP->solutions[condI]);

		// Fills row
		for (const auto& [condJ, bdrAttJ] : conductors) {
			// C_ij = Q_j / V_i. V_i is always 1.0
			C(condI, condJ) = sP->solver->getChargeInBoundary(conductors.at(condJ));
		}

		exportFieldSolutions(opts_, *sP->solver,
			model_.getMaterials().get<PEC>(condI).name, ignoreDielectrics);
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
	//   Clayton Paul's book: Analysis of Multiconductor Transmission Lines
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

void Driver::run()
{
	auto openness{ model_.determineOpenness() };
	if (openness == Model::Openness::closed) {
		auto pul = buildPULParametersForModel();
		saveToJSONFile(
			pul.toJSON(), 
			opts_.exportFolder + "pulmtln.out.json");
	}
	else if (openness == Model::Openness::open) {
		auto inCell = getInCellPotentials();
		saveToJSONFile(
			inCell.toJSON(),
			opts_.exportFolder + "inCellPotentials.out.json");
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

DenseMatrix Driver::getFloatingPotentialsMatrix(
	const bool ignoreDielectrics)
{
	// For an open-problem with N conductors, returns a NxN matrix which has: 
	// - a main diagonal of 1s, representing a prescribed voltage of 1 in the n-th conductor.
	// - the off-diagonal terms are the voltages at the other conductors when they are assumed to be floating.

	// For closed and semi-open problems with N conductors, returns a N-1 x N-1 matrix and assumes that conductor 0 has 
	// alway a prescribed potential of zero.

	if (model_.getMaterials().buildNameToAttrMapFor<PEC>().size() == 1) {
		mfem::DenseMatrix res(1,1);
		res = 1.0;
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
	mfem::DenseMatrix res(N, N);
	for (int i{ 0 }; i < N; ++i) {
		mfem::DenseMatrix A{ C };
		mfem::Vector negativeQ(N);
		negativeQ = 0.0;
		negativeQ(i) = -1.0;

		A.SetCol(i, negativeQ);

		mfem::Vector b(N);
		b = C.GetColumn(i);
		b *= -1.0;

		mfem::Vector x(N);
		mfem::DenseMatrixInverse Ainv(A);
		Ainv.Mult(b, x);

		for (int j{ 0 }; j < res.NumCols(); ++j) {
			if (i == j) {
				res(i, i) = 1.0;
			}
			else {
				res(i, j) = x(j);
			}
		}
	}

	return res;
}

FloatingPotentials Driver::getFloatingPotentials()
{
	FloatingPotentials res;

	res.electric = getFloatingPotentialsMatrix(false);
	res.magnetic = getFloatingPotentialsMatrix(true);

	return res;
}

std::list<std::string> listMaterialsInInnerRegion(
	const Model& m,
	bool includeConductors = true)
{
	std::list<std::string> res;
	res.push_back(INNER_REGION_DEFAULT_NAME);
	for (auto [name, tag] : m.getMaterials().buildNameToAttrMapFor<Dielectric>()) {
		if (name.find("Vacuum_") != std::string::npos) {
			continue;
		}
		res.push_back(name);
	}
	if (includeConductors) {
		for (auto [name, tag] : m.getMaterials().buildNameToAttrMapFor<PEC>()) {
			res.push_back(name);
		}
	}
	return res;
}


double Driver::getInnerRegionAveragePotential(
	const ElectrostaticSolver& s,
	bool includeConductors)
{

	double totalPotential = 0.0;
	double totalArea = 0.0;
	
	auto innerRegionMaterials = listMaterialsInInnerRegion(model_, includeConductors);
	auto materials = model_.getMaterials().buildNameToAttrMap();
	
	for (const auto& name: innerRegionMaterials ) {
		auto tag = materials.at(name);
		double area = model_.getAreaOfMaterial(name);
		if (model_.getMaterials().isDomainMaterial(name)) {
			totalPotential += s.getAveragePotentialInDomain(tag) * area;
		} else {
			totalPotential += s.getAveragePotentialInBoundary(tag) * area;
		}
		totalArea += area;
	}

	return totalPotential / totalArea;
}

std::map<MaterialId, FieldReconstruction> Driver::getFieldParameters(
	bool ignoreDielectrics)
{
	std::map<MaterialId, FieldReconstruction> res;

	auto fp = getFloatingPotentialsMatrix(ignoreDielectrics);

	SolvedProblem* sP;
	if (ignoreDielectrics) {
		sP = &magnetic_;
	}
	else {
		sP = &electric_;
	}

	ElectrostaticSolver& s = *sP->solver;

	const auto conductors{ model_.getMaterials().buildIdToAttrMapFor<PEC>() };
	for (const auto& [condI, bdrAttI] : conductors) {
		
		s.getPhi() *= 0.0;
		s.getE() *= 0.0;
		s.getD() *= 0.0;
		for (const auto& [condJ, bdrAttJ] : conductors) {
			s.getPhi().Add(fp(condI, condJ), *sP->solutions[condJ].phi);
			s.getE().Add(fp(condI, condJ), *sP->solutions[condJ].e);
			s.getD().Add(fp(condI, condJ), *sP->solutions[condJ].d);
		}

		exportFieldSolutions(opts_, s, 
			model_.getMaterials().get<PEC>(condI).name + "_prescribed_and_others_floating", 
			ignoreDielectrics);

		res[condI].innerRegionAveragePotential = 
			getInnerRegionAveragePotential(s, true);
		auto centerOfCharge = s.getCenterOfCharge();
		std::copy(
			centerOfCharge.begin(), centerOfCharge.end(), 
			res[condI].expansionCenter.begin());
		res[condI].ab = s.getMultipolarCoefficients(opts_.multipolarExpansionOrder);
		for (const auto& [condJ, bdrAttJ] : conductors) {
			res[condI].conductorPotentials[condJ] = fp(condI, condJ);
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

	res.innerRegionBox = model_.getBoundingBoxOfMaterial(INNER_REGION_DEFAULT_NAME);

	res.electric = getFieldParameters(false);
	res.magnetic = getFieldParameters(true);

	return res;
}

}