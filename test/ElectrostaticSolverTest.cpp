#include <gtest/gtest.h>

#include <functional>

#include "TestUtils.h"
#include "ElectrostaticSolver.h"
#include "Parser.h"
#include "Driver.h"

using namespace mfem;
using namespace pulmtln;

Vector getBaricenterOfElement(Mesh& mesh, int e)
{
	Element* elem{ mesh.GetElement(e) };
	mfem::Vector center({ 0.0, 0.0 });
	mfem::Array<int> vs(elem->GetNVertices());
	elem->GetVertices(vs);
	for (auto i{ 0 }; i < vs.Size(); i++) {
		mfem::Vector vertexPos(mesh.GetVertex(vs[i]), 2);
		center += vertexPos;
	}
	center /= vs.Size();
	return center;
}

void wireField(const Vector& pos, Vector& res)
{
	// Assumes a charge of 1 unit in the wire.
	double norm = pos.Norml2();
	res = pos;
	res /= (norm * 2.0 * M_PI);
}

void exportSolution(ElectrostaticSolver& s, const std::string& caseName)
{
	ParaViewDataCollection paraview_dc{ outFolder() + caseName, s.getMesh() };
	s.writeParaViewFields(paraview_dc);
}

class ElectrostaticSolverTest : public ::testing::Test {
};


TEST_F(ElectrostaticSolverTest, parallel_plates)
{
	//    0 V
	//  +-----+
	//  |  3  |
	//  |4   2|
	//  |  1  |
	//  +-----+
	//    1 V
	auto m{ Mesh::MakeCartesian2D(1, 5, Element::QUADRILATERAL, 1.0, 1.0) };

	SolverInputs params;
	params.dirichletBoundaries = {
		{
			{1,    1.0}, // bottom boundary.
			{3,    0.0}, // top boundary.
		} 
	};

	ElectrostaticSolver s{m, params};
	s.Solve();

	exportSolution(s, "parallel_plates");

	const double aTol{ 1e-5 };
	const double rTol{ 1e-5 }; // Error percentage of 0.001%

	EXPECT_NEAR(0.0, s.getTotalCharge(), aTol);

	EXPECT_LE(relError(1.0, s.getChargeInBoundary(1)), rTol);

	EXPECT_LE(relError(-1.0, s.getChargeInBoundary(3)), rTol);

	EXPECT_NEAR(0.0, s.getChargeInBoundary(2), aTol);

	// Expected energy formula (1/2) * epsilon_0 * E^2 * A.
	// Area A = 1.0 * 1.0
	// Electric field is constant = 1.0 V/m.
	double expectedEnergy{ 0.5 * EPSILON0_NATURAL };
	EXPECT_LE(relError(expectedEnergy, s.getTotalEnergy()), rTol);

}

TEST_F(ElectrostaticSolverTest, parallel_plates_energy)
{
	//    0 V
	//  +-----+
	//  |  3  |
	//  |4   2|
	//  |  1  |
	//  +-----+
	//    1 V
	auto m{ Mesh::MakeCartesian2D(1, 5, Element::QUADRILATERAL, 1.0, 1.0) };

	SolverInputs params;
	params.dirichletBoundaries = {
		{
			{1,    1.0}, // bottom boundary.
			{3,    0.0}, // top boundary.
		}
	};

	ElectrostaticSolver s{ m, params };
	s.Solve();

	// Expected energy formula (1/2) * epsilon_0 * E^2 * A.
	// Area A = 1.0 * 1.0
	// Electric field is constant = 1.0 V/m.
	const double rTol{ 1e-5 }; // Error percentage of 0.001%
	double expectedEnergy{ 0.5 * EPSILON0_NATURAL };
	EXPECT_LE(relError(expectedEnergy, s.getTotalEnergy()), rTol);

}

TEST_F(ElectrostaticSolverTest, parallel_plates_neumann)
{
	//    Q = -1 
	//  +-----+
	//  |  3  |
	//  |4   2|
	//  |  1  |
	//  +-----+
	//    Q = 1 
	auto m{ Mesh::MakeCartesian2D(5, 5, Element::QUADRILATERAL, 1.0, 1.0) };

	SolverInputs params;
	params.neumannBoundaries = {
		{
			{1,  1.0}, 
			{3, -1.0}, 
		}
	};

	ElectrostaticSolver s{ m, params };
	s.Solve();

	exportSolution(s, "parallel_plates_neumann");


	const double rTol{ 5e-2 }; 
	EXPECT_LE(relError( 1.0, s.getChargeInBoundary(1)), rTol);
    EXPECT_LE(relError(-1.0, s.getChargeInBoundary(3)), rTol);

	const double aTol{ 1e-5 };
	EXPECT_NEAR(0.0, s.getTotalCharge(), aTol);
	EXPECT_NEAR(0.0, s.getChargeInBoundary(2), aTol);
	EXPECT_NEAR(0.0, s.getChargeInBoundary(4), aTol);

}

TEST_F(ElectrostaticSolverTest, parallel_plates_epsr2)
{
	//    0 V
	//  +-----+
	//  |  3  |
	//  |4   2|
	//  |  1  |
	//  +-----+
	//    1 V
	auto m{ Mesh::MakeCartesian2D(1, 5, Element::QUADRILATERAL, 1.0, 1.0) };

	SolverInputs p;
	p.dirichletBoundaries = { {
		{1,    1.0}, // bottom boundary.
		{3,    0.0}, // top boundary.
	} };
	p.domainPermittivities = {{ {1, 2.0} }};

	ElectrostaticSolver s{m, p, SolverOptions{}};
	s.Solve();

	exportSolution(s, "Parallel_plates_epsr2");

	const double tol{ 1e-6 };
	EXPECT_LE(relError(2.0, s.getChargeInBoundary(1)), tol);
}

TEST_F(ElectrostaticSolverTest, two_materials)
{
	//         0 V,   Bdr 3
	//       +-------------+
	//       | eps_r=1, D1 |
	// Bdr 4 |-------------| Bdr 2
	//       | eps_r=4, D2 |
	//       +-------------+
	//         1 V,   Bdr 1
	auto m{ Mesh::MakeCartesian2D(1, 10, Element::QUADRILATERAL, 1.0, 1.0) };
	for (auto i{ 0 }; i < m.GetNE(); ++i) {
		auto center{ getBaricenterOfElement(m, i) };
		if (center[1] < 0.5) {
			m.SetAttribute(i, 2);
		}
	}
	m.Finalize();
	m.SetAttributes();

	SolverInputs p;
	p.dirichletBoundaries = { {
		{1, 1.0}, // bottom boundary.
		{3, 0.0}, // top boundary.
	} };
	p.domainPermittivities = { {{2, 4.0}} };

	ElectrostaticSolver s{m, p};
	s.Solve();

	exportSolution(s, "two_materials");

	const double aTol{ 1e-4 };
	const double rTol{ 1e-5 };  // 0.001%

	EXPECT_NEAR(0.0, s.getTotalCharge(), aTol);

	EXPECT_LE(relError(1.6, s.getChargeInBoundary(1)), rTol);
	EXPECT_LE(relError(- 1.6, s.getChargeInBoundary(3)), rTol);

}

TEST_F(ElectrostaticSolverTest, empty_coax_charge_in_boundaries)
{
	// Coaxial case.
	const std::string CASE{ "empty_coax" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	const double V{ 1.0 };
	SolverInputs p;
	p.dirichletBoundaries = {{
		{1, 0.0}, // outer boundary
		{2, V},   // inner boundary
	}};
	
	ElectrostaticSolver s{m, p};
	s.Solve();

	exportSolution(s, getCaseName());

	// Expected capacitance C = eps0 * 2 * pi / log(ro/ri)
	// Expected charge      QExpected = V0 * C 
	double QExpected{ V * EPSILON0_NATURAL * 2 * M_PI / log(0.05 / 0.025) };

	const double rTol{ 5e-3 }; // 0.5% error.

	EXPECT_LE(relError( QExpected, s.getChargeInBoundary(2)), rTol); // Boundary 2 is the internal.
	EXPECT_LE(relError(-QExpected, s.getChargeInBoundary(1)), rTol); // Boundary 1 is the external.
}

TEST_F(ElectrostaticSolverTest, empty_coax_average_potential)
{
	// Coaxial case.
	const std::string CASE{ "empty_coax" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	SolverInputs p;
	p.dirichletBoundaries = { {
		{2, 1.0},   // inner boundary
	} };
	p.openBoundaries = { 1 };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	exportSolution(s, getCaseName());

	double avPotential = s.getAveragePotentialInDomain(3);

	const double a = 0.05;  // external boundary radius
	const double b = 0.025; // inner boundary radius
	double Qi = s.getChargeInBoundary(2);
	double area = M_PI * (a * a - b * b);
	double totalVoltage =
		( 1.0 + Qi/2.0/M_PI/EPSILON0_NATURAL*std::log(b) ) * area
		- (Qi / 2.0 / EPSILON0_NATURAL) * (a * a * (std::log(a)-0.5) - b * b * (std::log(b)-0.5));
		
	double expected = totalVoltage / area;

	EXPECT_NEAR(0.0, relError(expected, avPotential), 5e-4);
}


TEST_F(ElectrostaticSolverTest, empty_coax_neumann)
{
	// Coaxial case.
	const std::string CASE{ "empty_coax" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	const double V{ 1.0 };
	SolverInputs p;
	p.neumannBoundaries = { {
		{2, 1.0 / (2 * M_PI * 25e-3)},   // Inner boundary, 1 unit of charge in total.
	} };
	p.dirichletBoundaries = { {
		{1, 0.0},   // Outer boundary. 
	} };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	exportSolution(s, getCaseName());

	const double rTol{ 5e-3 }; // 0.5% error.
	EXPECT_LE(relError( 1.0, s.getChargeInBoundary(2)), rTol); // Boundary 2 is the internal.
	EXPECT_LE(relError(-1.0, s.getChargeInBoundary(1)), rTol); // Boundary 1 is the external.
}

TEST_F(ElectrostaticSolverTest, empty_coax_neumann_quadrupole)
{
	// Coaxial case.
	const std::string CASE{ "empty_coax" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	SolverInputs p;
	
	std::vector<multipolarCoefficient> ab = {
		{0.0, 0.0}, // a0
		{0.0, 0.0}, // a1, b1
		{1.0, 0.0}, // a2, b2
	};
	Vector origin({ 0.0, 0.0 });
	
	ElectrostaticSolver s{ m, p };

	// Sets multipolar expansion over internal boundary.
	{
		std::function<double(const Vector&)> f =
			std::bind(&multipolarExpansion, std::placeholders::_1, ab, origin);
		FunctionCoefficient fc(f);
		s.setNeumannCondition(2, fc);
	}

	// Sets opposite multipolar expansion over external boundary.
	{
		ab[2] = { -1.0, 0.0 };
		std::function<double(const Vector&)> f =
			std::bind(&multipolarExpansion, std::placeholders::_1, ab, origin);
		FunctionCoefficient fc(f);
		s.setNeumannCondition(1, fc);
	}

			
	s.Solve();

	exportSolution(s, getTestCaseName());

	const double aTol{ 5e-4 }; 
	EXPECT_NEAR(0.0, s.getChargeInBoundary(2), aTol); // Boundary 2 is the internal.
	EXPECT_NEAR(0.0, s.getChargeInBoundary(1), aTol); // Boundary 1 is the external.
}

TEST_F(ElectrostaticSolverTest, wire_in_open_region)
{
	// Coaxial case.
	const std::string CASE{ "empty_coax" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };
	
	SolverInputs p;
	p.dirichletBoundaries = { {
		{2, 1.0}
	} };
	p.openBoundaries = { 1 }; // Outer boundary.

	ElectrostaticSolver s{ m, p };
	s.Solve();

	exportSolution(s, getCaseName());

	double Q = s.getChargeInBoundary(2);
	ConstantCoefficient chargeCoeff{ Q };
	VectorFunctionCoefficient exactField{2, wireField, &chargeCoeff};
	auto error = s.getE().ComputeL2Error(exactField);
	EXPECT_LE(error, 1.2);

	double CExpected = EPSILON0_NATURAL * 2 * M_PI / log(0.05 / 0.025);
	const double rTol{ 1e-2 }; // 1% error.

	auto U{ s.getTotalEnergy()};
	double CFromEnergy = 0.5 * std::pow(Q,2) / U;
	EXPECT_LE(relError(CExpected, CFromEnergy), rTol);

	double Qb = s.getChargeInBoundary(1);
	EXPECT_LE(relError(Q, -Qb), rTol);

	double Vb = s.getAveragePotentialInBoundary(1);
	auto CFromVb = EPSILON0_NATURAL * Q / (1.0 - Vb);
	EXPECT_LE(relError(CExpected, CFromVb), rTol);

	exportSolution(s, getCaseName());

}

TEST_F(ElectrostaticSolverTest, two_wires_coax)
{
	const std::string CASE{ "two_wires_coax" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	const double V{ 1.0 }; // Voltage
	SolverInputs p;
	p.dirichletBoundaries = { {
		{1, 0.0}, // Conductor 0 bdr (GND).
		{2, V},   // Conductor 1 bdr.
		{3, 0.0}, // Conductor 2 bdr.
	} };

	ElectrostaticSolver s{m, p};
	s.Solve();

	exportSolution(s, getCaseName());

	double C11Expected{  2.15605359 };
	double C12Expected{ -0.16413431 };

	const double rTol{ 2.5e-2 };

	EXPECT_LE(relError(C11Expected, s.getChargeInBoundary(2) / V), rTol);
	EXPECT_LE(relError(C12Expected, s.getChargeInBoundary(3) / V), rTol);
}

TEST_F(ElectrostaticSolverTest, two_wires_open_capacitance)
{
	const std::string CASE{ "two_wires_open" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	const double V{ 1.0 }; // Voltage
	SolverInputs p;
	p.dirichletBoundaries = { {
		{1,  V}, // Conductor 1 bdr.
		{2, -V}, // Conductor 2 bdr.
	} };
	p.openBoundaries = { 3 };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	exportSolution(s, getCaseName());

	// Capacity between two straight wires.
	double d = 50;
	double rw1 = 2;
	double rw2 = 2;

	double CExpected{ 
		2*M_PI*EPSILON0_NATURAL /
		std::acosh( (d*d - rw1*rw1 - rw2*rw2) / (2*rw1*rw2) )
	};
	
	double chargeInOpenBoundary{ s.getChargeInBoundary(3) };
	EXPECT_LE(1e-6, std::abs(chargeInOpenBoundary));

	const double rTol{ 5e-4 };
	double CComputed{ s.getChargeInBoundary(1) / (2*V) };
	EXPECT_LE(relError(CExpected, CComputed), rTol);

	double CComputedEnergy{ 2.0 * s.getTotalEnergy() / (4.0*V*V) };
	EXPECT_LE(relError(CExpected, CComputedEnergy), rTol);
}

TEST_F(ElectrostaticSolverTest, two_wires_open_monopolar_moment)
{
	const std::string CASE{ "two_wires_open" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	SolverInputs p;
	p.openBoundaries = { 3 };
	p.dirichletBoundaries = { {
		{1,  1.0}, // Conductor 1 bdr.
		{2,  1.0}, // Conductor 2 bdr.
	} };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	auto Q1{ s.getChargeInBoundary(1) };
	auto Q2{ s.getChargeInBoundary(2) };
	auto Qt{ Q1 + Q2 };

	auto a0 = s.getChargeMomentComponent(0, 0, Vector({ 0.0, 0.0 }));
	EXPECT_NEAR(Qt, a0, 5e-6);

	auto b0 = s.getChargeMomentComponent(0, 1, Vector({ 0.0, 0.0 }));
	EXPECT_NEAR(0.0, b0, 5e-6);

	exportSolution(s, getCaseName());
}

TEST_F(ElectrostaticSolverTest, two_wires_open_boundary_charges)
{
	const std::string CASE{ "two_wires_open" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	const double V{ 1.0 }; // Voltage
	SolverInputs p;
	p.dirichletBoundaries = { {
		{1,  V}, // Conductor 1 bdr.
		{2,  V}, // Conductor 2 bdr.
	} };
	p.openBoundaries = { 3 };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	exportSolution(s, getTestCaseName());


	auto Q1{ s.getChargeInBoundary(1) };
	auto Q2{ s.getChargeInBoundary(2) };
	auto Qb{ s.getChargeInBoundary(3) };

	auto Vb{ s.getAveragePotentialInBoundary(3) };
	auto Vd = V - Vb;

	EXPECT_NEAR(0.0, Q1 + Q2 + Qb, 1e-3);
	
	auto CFromEnergy{ 2.0 * s.getTotalEnergy() / (Vd * Vd) };
	
	auto C1b = Q1 / Vd;
	auto C2b = Q2 / Vd;
	auto CFromCharge{ C1b + C2b };

	EXPECT_LE(relError(CFromEnergy, CFromCharge), 1e-3);
}

TEST_F(ElectrostaticSolverTest, two_wires_open_multipolarCoefficients_with_same_potential)
{
	const std::string CASE{ "two_wires_open" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	SolverInputs p;
	p.openBoundaries = { 3 };
	
	p.dirichletBoundaries = { {
		{1,  1.0}, // Conductor 1 bdr.
		{2,  1.0}, // Conductor 2 bdr.
	} };

	ElectrostaticSolver s{ m, p };
	s.Solve();
		
	auto Q1 = s.getChargeInBoundary(1);
	auto Q2 = s.getChargeInBoundary(2);
	auto Qt = Q1 + Q2;

	auto ab = s.getMultipolarCoefficients(2); // Up to quadrupolar moment.
	
	ASSERT_EQ(3, ab.size());

	const double aTol = 1e-5;
	EXPECT_NEAR(Qt, ab[0].first,  aTol);
	EXPECT_NEAR(0.0, ab[0].second, aTol);

	EXPECT_NEAR(0.0, ab[1].first, aTol);
	EXPECT_NEAR(0.0, ab[1].second, aTol);

	const double a = 0.025;
	const double quadrupoleTerm = std::pow(a, 2) * Q1;
	EXPECT_NEAR(quadrupoleTerm, ab[2].first, aTol);
	EXPECT_NEAR(0.0, ab[2].second, aTol);
}

TEST_F(ElectrostaticSolverTest, two_wires_open_center_of_charge_with_same_potential)
{
	const std::string CASE{ "two_wires_open" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	SolverInputs p;
	p.openBoundaries = { 3 };
	p.dirichletBoundaries = { {
		{1,  1.0}, // Conductor 1 bdr.
		{2,  1.0}, // Conductor 2 bdr.
	} };

	ElectrostaticSolver s{ m, p };
	s.Solve();
	auto centerOfCharge{ s.getCenterOfCharge() };
	EXPECT_NEAR(0.0, centerOfCharge[0], 5e-6);
	EXPECT_NEAR(0.0, centerOfCharge[1], 5e-6);

	exportSolution(s, getCaseName());
}

TEST_F(ElectrostaticSolverTest, two_wires_open_center_of_charge_with_floating_potential)
{
	const std::string CASE{ "two_wires_open" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	SolverInputs p;
	p.openBoundaries = { 3 };

	p.dirichletBoundaries = { {
		{1,  1.0}, // Conductor 1 bdr.
		{2,  0.48228164}, // Conductor 2, floating, zero charge.
	} };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	ASSERT_NEAR(0.0, s.getChargeInBoundary(2), 5e-5);

	auto centerOfCharge{ s.getCenterOfCharge() };
	EXPECT_NEAR(-0.025, centerOfCharge[0], 1e-4);
	EXPECT_NEAR(0.0, centerOfCharge[1], 1e-4);

	exportSolution(s, getCaseName());
}

TEST_F(ElectrostaticSolverTest, three_wires_ribbon_zero_net_charge)
{
	// Three wires ribbon open problem. 
	// Comparison with Clayton Paul's book:  
	// Analysis of multiconductor transmision lines. 2007.
	// Sec. 5.2.3, p. 187.

	const std::string CASE{ "three_wires_ribbon" };
	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	SolverInputs p;
	p.dirichletBoundaries = {
		{
			{1, 1.0}, // Conductor 0
			{2, 0.0}, // Conductor 1
			{3, 1.0}, // Conductor 2
		}
	};

	SolverOptions solverOpts;
	ElectrostaticSolver s{ m, p, solverOpts };
	s.Solve();

	exportSolution(s, getCaseName());

	auto Q0 = s.getChargeInBoundary(1);
	auto Q1 = s.getChargeInBoundary(2);
	auto Q2 = s.getChargeInBoundary(3);
	auto Qb = s.getChargeInBoundary(4);
	EXPECT_NEAR(0.0, Q0+Q1+Q2+Qb, 1e-4);
}

TEST_F(ElectrostaticSolverTest, lansink2024_fdtd_in_cell_C00_with_floating)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.

	const std::string CASE{ "lansink2024_fdtd_cell" };
	const std::string fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
	auto model{ Parser{fn}.readModel() };

	auto fp = Driver::loadFromFile(fn).getFloatingPotentials().electric;

	SolverInputs p;
	p.dirichletBoundaries = {
		{
			{1, fp(0,0)}, // Conductor 0,
			{2, fp(0,1)}  // Conductor 1
		}
	};
	p.openBoundaries = { 3 };

	ElectrostaticSolver s{ *model.getMesh(), p };
	s.Solve();

	auto avVVacuum = s.getAveragePotentialInDomain(5);
	auto avV0 = s.getAveragePotentialInBoundary(1);
	auto avV1 = s.getAveragePotentialInBoundary(2);

	auto Q0 = s.getChargeInBoundary(1);
	auto Q1 = s.getChargeInBoundary(2);

	auto areaVacuum = model.getAreaOfMaterial("Vacuum_0");
	auto areaCond0 = model.getAreaOfMaterial("Conductor_0");
	auto areaCond1 = model.getAreaOfMaterial("Conductor_1");
	auto totalArea = areaVacuum + areaCond0 + areaCond1;

	auto avV =
		(avVVacuum * areaVacuum + avV0 * areaCond0 + avV1 * areaCond1) / (totalArea);
	avV = -avV + fp(0, 0); // In the paper, Conductor_0 is assumed to have zero voltage.

	auto computedC00 = std::abs(Q0 / avV * EPSILON0_SI); // C11 in paper.

	exportSolution(s, getTestCaseName());

	// from Table 1, floating conductor case. C11
	auto expectedC00 = 14.08e-12;

	// 
	double rTol = 0.01;
	EXPECT_NEAR(0.0, relError(expectedC00, computedC00), rTol);
	EXPECT_NEAR(0.0, Q1, 0.005);
}

TEST_F(ElectrostaticSolverTest, lansink2024_fdtd_in_cell_C01_with_floating)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.

	const std::string CASE{ "lansink2024_fdtd_cell" };
	const std::string fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
	auto model{ Parser{fn}.readModel() };

	auto fp = Driver::loadFromFile(fn).getFloatingPotentials().electric;

	SolverInputs p;
	p.dirichletBoundaries = {
		{
			{1, fp(1,0)}, // Conductor 0,
			{2, fp(1,1)}  // Conductor 1
		}
	};
	p.openBoundaries = { 3 };

	ElectrostaticSolver s{ *model.getMesh(), p };
	s.Solve();

	auto avVVacuum = s.getAveragePotentialInDomain(5);
	auto avV0 = s.getAveragePotentialInBoundary(1);
	auto avV1 = s.getAveragePotentialInBoundary(2);

	auto Q0 = s.getChargeInBoundary(1);
	auto Q1 = s.getChargeInBoundary(2);

	auto areaVacuum = model.getAreaOfMaterial("Vacuum_0");
	auto areaCond0 = model.getAreaOfMaterial("Conductor_0");
	auto areaCond1 = model.getAreaOfMaterial("Conductor_1");
	auto totalArea = areaVacuum + areaCond0;

	auto avV =
		(avVVacuum * areaVacuum + avV0 * areaCond0) / (totalArea);
	avV = -avV + fp(1,0); // In the paper, Conductor_0 is assumed to have zero voltage.

	auto computedC01 = std::abs(Q1 / avV * EPSILON0_SI);

	exportSolution(s, getTestCaseName());

	// from Table 1, floating conductor case. C12 in paper
	auto expectedC01 = 43.99e-12; 

	// 
	double rTol = 0.025;
	EXPECT_NEAR(0.0, relError(expectedC01, computedC01), rTol);
	EXPECT_NEAR(0.0, Q0, 0.005);
}

TEST_F(ElectrostaticSolverTest, lansink2024_single_wire_L00_with_floating)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.

	const std::string CASE{ "lansink2024_single_wire" };
	const std::string fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
	auto model{ Parser{fn}.readModel() };

	auto fp = Driver::loadFromFile(fn).getFloatingPotentials().electric;

	SolverInputs p;
	p.dirichletBoundaries = {
		{
			{1, 1.0}, // Conductor 0,
		}
	};
	p.openBoundaries = { 2 };

	ElectrostaticSolver s{ *model.getMesh(), p };
	s.Solve();

	auto avVVacuum = s.getAveragePotentialInDomain(4);
	auto avVDielectric = s.getAveragePotentialInDomain(5);
	auto avV0 = s.getAveragePotentialInBoundary(1);

	auto Q0 = s.getChargeInBoundary(1);

	auto areaVacuum = model.getAreaOfMaterial("Vacuum_0");
	auto areaDielectric = model.getAreaOfMaterial("Dielectric_0");
	auto areaCond0 = model.getAreaOfMaterial("Conductor_0");
	auto totalArea = areaVacuum + areaDielectric + areaCond0;

	auto avV =
		(avVVacuum * areaVacuum + avV0 * areaCond0 + avVDielectric * areaDielectric) / (totalArea);
	avV = -avV + 1.0;
	
	auto computedL00 = avV / Q0 * MU0_SI; // L11 in paper 

	exportSolution(s, getTestCaseName());

	// Table 3 result has a mistake. 
	// This is the correct value obtained through personal communication.
	auto expectedL00 = 320e-9;

	// 
	double rTol = 0.04;
	EXPECT_NEAR(0.0, relError(expectedL00, computedL00), rTol);
}

TEST_F(ElectrostaticSolverTest, lansink2024_small_one_centered_bem_comparison)
{
	// Case from:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.

	// SMALL WIRE IS CENTERED AT ORIGIN.
	// This test compares multipolar expansion results between Tulip and BEM.

	const std::string CASE{ "lansink2024_small_one_centered" };
	const std::string fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
	auto model{ Parser{fn}.readModel() };

	auto fp = Driver::loadFromFile(fn).getFloatingPotentials().electric;

	SolverInputs p;
	p.dirichletBoundaries = {
		{
			{1, fp(0,0)}, // Conductor 0,
			{2, fp(0,1)}, // Conductor 1,
		}
	};
	p.openBoundaries = { 3 };

	ElectrostaticSolver s{ *model.getMesh(), p };
	s.Solve();
	

	int order = 2;

	// Tulip coefficients.
	multipolarCoefficients ab(order + 1);
	mfem::Vector origin({ 0.0, 0.0 });
	for (int n = 0; n < order + 1; n++) {
		ab[n] = {
			s.getChargeMomentComponent(n, 0, origin),
			s.getChargeMomentComponent(n, 1, origin)
		};
	}

	// BEM coefficients.
	multipolarCoefficients abBEM(order + 1);
	abBEM[0] = { 1.0, 0.0 };
	abBEM[1] = { 0.005, 0.0 };
	abBEM[2] = { 8.7819e-5, 0 };

	// Checks normalizing by total charge.
	const double aTol = 1e-4;
	auto a0 = ab[0].first;
	for (int n = 0; n < order + 1; n++) {
		EXPECT_NEAR(ab[n].first / a0,  abBEM[n].first,  aTol) << "for n=" << n;
		EXPECT_NEAR(ab[n].second / a0, abBEM[n].second, aTol) << "for n=" << n;
	}

	// For debugging.
	// Export FEM solution.
	exportSolution(s, getTestCaseName());

	// Export multipolar expansion projection.
	auto multipolarSolution = s.getSolution();
	std::function<double(const Vector&)> f =
		std::bind(&multipolarExpansion, std::placeholders::_1, ab, origin);
	FunctionCoefficient fc(f);
	multipolarSolution.phi->ProjectCoefficient(fc);
	s.setSolution(multipolarSolution);
	exportSolution(s, getTestCaseName()+"_from_multipolar_expansion_in_origin");

	// Computes and exports multipolar expansion with center of charge.
	{
		ElectrostaticSolver s{ *model.getMesh(), p };
		s.Solve();
		auto centerOfCharge = s.getCenterOfCharge();
		auto mCoeff = s.getMultipolarCoefficients(order);
		std::function<double(const Vector&)> f =
			std::bind(&multipolarExpansion, std::placeholders::_1, mCoeff, centerOfCharge);
		FunctionCoefficient fc(f);
		auto multipolarSolution = s.getSolution();
		multipolarSolution.phi->ProjectCoefficient(fc);
		s.setSolution(multipolarSolution);
		exportSolution(s, getTestCaseName() + "_from_multipolar_expansion_in_center_of_charge");
	}

}