#include "gtest/gtest.h"

#include "TestUtils.h"

#include "constants.h"
#include "Driver.h"
#include "Solver.h"

using namespace tulip;

using json = nlohmann::json;

class DriverTest : public ::testing::Test {};

TEST_F(DriverTest, empty_coax)
{
	// Empty Coaxial case.
	auto out{ Driver::loadFromAdaptedFile(inputCase("empty_coax")).getPULMTL() };

	auto CExpected{ EPSILON0_SI * 2 * M_PI / log(0.05 / 0.025) };

	const double rTol{ 0.005 };
	ASSERT_EQ(1, out.C.NumCols() * out.C.NumRows());
	EXPECT_LE(relError(CExpected, out.C(0, 0)), rTol);

	auto LExpected{ EPSILON0_SI * MU0_SI / CExpected };
	ASSERT_EQ(1, out.L.NumCols() * out.L.NumRows());
	EXPECT_LE(relError(LExpected, out.L(0, 0)), rTol);
}

TEST_F(DriverTest, partially_filled_coax)
{
	// Partially filled coax.
	// External radius -> r0 = 50 mm
	// Internal radius -> r1 = 25 mm 
	// Dielectric internal radius -> rI_dielectric = 25 mm
	// Dielectric external radius -> rO_dielectric = 35 mm
	// Dielectric permittivity -> eps_r = 4.0
	auto out{ Driver::loadFromAdaptedFile(inputCase("partially_filled_coax")).getPULMTL() };

	// Equivalent capacity is the series of the inner and outer capacitors.
	auto COut{ EPSILON0_SI * 2 * M_PI / log(0.050 / 0.035) };
	auto CIn{ 4.0 * EPSILON0_SI * 2 * M_PI / log(0.035 / 0.025) };
	auto CExpected = COut * CIn / (COut + CIn);

	const double rTol{ 0.002 }; // 0.2% Error.
	ASSERT_EQ(1, out.C.NumCols() * out.C.NumRows());
	EXPECT_LE(relError(CExpected, out.C(0, 0)), rTol);

	// Inductance can be calculated from free-space capacity (C0).
	auto C0 = EPSILON0_SI * 2 * M_PI / log(0.05 / 0.025);
	auto LExpected{ EPSILON0_SI * MU0_SI / C0 };

	ASSERT_EQ(1, out.L.NumCols() * out.L.NumRows());
	EXPECT_LE(relError(LExpected, out.L(0, 0)), rTol);
}

TEST_F(DriverTest, partially_filled_coax_by_domains)
{
	auto dr{ Driver::loadFromAdaptedFile(inputCase("partially_filled_coax")) };
	auto globalOut{ dr.getPULMTL() };

	auto out{ dr.getPULMTLByDomains() };
	EXPECT_EQ(1, out.domainTree.verticesSize());
	ASSERT_EQ(1, out.domainToPUL.size());


	// There are some minor differences in output. 
	// I do not know why but my guess is that it is due to initial seeds in
	// the iterative solver.
	const double rTol{ 1e-6 };
	ASSERT_EQ(1, out.domainToPUL.count(0));
	ASSERT_EQ(globalOut.L.NumRows(), out.domainToPUL.at(0).L.NumRows());
	ASSERT_EQ(globalOut.C.NumRows(), out.domainToPUL.at(0).C.NumRows());
	EXPECT_LE(relError(globalOut.L(0, 0), out.domainToPUL.at(0).L(0, 0)), rTol);
	EXPECT_LE(relError(globalOut.C(0, 0), out.domainToPUL.at(0).C(0, 0)), rTol);
}

TEST_F(DriverTest, two_wires_coax)
{
	const std::string CASE{ "two_wires_coax" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };

	mfem::DenseMatrix CMatExpected(2, 2);
	CMatExpected(0, 0) = 2.15605359;
	CMatExpected(0, 1) = -0.16413431;
	CMatExpected(1, 0) = CMatExpected(0, 1);
	CMatExpected(1, 1) = CMatExpected(0, 0);
	CMatExpected *= EPSILON0_SI;

	const double rTol{ 2.5e-2 };

	auto out{ Driver::loadFromAdaptedFile(fn).getPULMTL() };

	const int N{ 2 };
	ASSERT_EQ(N, out.C.NumCols());
	ASSERT_EQ(N, out.C.NumRows());

	// Compares with analytical solution.
	for (int i{ 0 }; i < N; i++) {
		for (int j{ 0 }; j < N; j++) {
			EXPECT_LE(relError(CMatExpected(i, j), out.C(i, j)), rTol);
		}
	}

	// Checks matrix are symmetric.
	for (int i{ 0 }; i < N; i++) {
		for (int j{ 0 }; j < N; j++) {
			EXPECT_EQ(out.C(i, j), out.C(j, i));
			EXPECT_EQ(out.L(i, j), out.L(j, i));
		}
	}
}

TEST_F(DriverTest, two_wires_shielded_floating_potentials)
{
	const std::string CASE{ "two_wires_shielded" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };

	auto driver = Driver::loadFromAdaptedFile(fn); 

	ASSERT_ANY_THROW(driver.getFloatingPotentials(1,false));
}

TEST_F(DriverTest, two_wires_open_floating_potentials)
{
	const std::string CASE{ "two_wires_open" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };

	auto driver{ Driver::loadFromAdaptedFile(fn) };
	auto fp{ driver.getFloatingPotentials(0, false) };

	ASSERT_EQ(2, fp.size());

	// Verifies that floating potentials result in zero charge on the floating conductor.
	auto eSol = driver.getElectrostaticSolvedProblem();
	driver.loadFloatingPotentials(eSol, fp);

	// For debugging.
	auto& s = *eSol->solver;
	ParaViewDataCollection pd(outFolder() + CASE + "_floating", s.getMesh());
	s.writeParaViewFields(pd);
	
	auto cond1 = driver.getModel().getMaterials().getConductorWithId(1);
	auto Q1 = s.getChargeInBoundary(cond1->getAttribute());

	EXPECT_NEAR(0.0, Q1, 1e-4); // Floating conductor should have zero charge.
}

TEST_F(DriverTest, five_wires)
{
	// Five wires in round shield. 
	// Comparison with SACAMOS data (No Laplace).

	const std::string CASE{ "five_wires" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };

	mfem::DenseMatrix couplingExpected(5, 5);
	double couplingExpectedData[25] = {
		 1.0000, -0.2574    , -0.2574    , -0.2574    , -0.2574    ,
		-0.2574,  1.0000    , -0.029017  , -3.5871E-05, -0.029017  ,
		-0.2574, -0.029017  ,  1         , -0.029017  , -3.5871E-05,
		-0.2574, -3.5871E-05, -0.029017	 ,  1         , -0.029017  ,
		-0.2574, -0.029017  , -3.5871E-05, -0.029017  ,  1.0000
	};
	couplingExpected.UseExternalData(couplingExpectedData, 5, 5);

	auto out{
		Driver::loadFromAdaptedFile(fn).getPULMTL().getCapacitiveCouplingCoefficients()
	};

	ASSERT_EQ(couplingExpected.NumRows(), out.NumRows());
	ASSERT_EQ(couplingExpected.NumCols(), out.NumCols());

	double rTol{ 0.05 };
	for (int i{ 0 }; i < couplingExpected.NumRows(); i++) {
		for (int j{ 0 }; j < couplingExpected.NumCols(); j++) {
			EXPECT_LE(std::abs(couplingExpected(i, j) - out(i, j)), rTol)
				<< "In C(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, three_wires_ribbon)
{
	// Three wires ribbon open problem. 
	// Comparison with Clayton Paul's book:  
	// Analysis of multiconductor transmision lines. 2007.
	// Sec. 5.2.3, p. 187.

	const std::string CASE{ "three_wires_ribbon" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };

	double CExpectedData[4] = {
		 37.432, -18.716,
		 -18.716, 24.982
	};
	mfem::DenseMatrix CExpected(2, 2);
	CExpected.UseExternalData(CExpectedData, 2, 2);
	CExpected *= 1e-12;

	double LExpectedData[4] = {
		0.74850, 0.50770,
		0.50770, 1.0154
	};
	mfem::DenseMatrix LExpected(2, 2);
	LExpected.UseExternalData(LExpectedData, 2, 2);
	LExpected *= 1e-6;

	auto out{ Driver::loadFromAdaptedFile(fn).getPULMTL() };

	const double rTol{ 0.005 };
	ASSERT_EQ(CExpected.NumRows(), out.C.NumRows());
	ASSERT_EQ(CExpected.NumCols(), out.C.NumCols());
	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
			EXPECT_LE(relError(CExpected(i, j), out.C(i, j)), rTol) <<
				"In C(" << i << ", " << j << ")";
		}
	}

	ASSERT_EQ(LExpected.NumRows(), out.L.NumRows());
	ASSERT_EQ(LExpected.NumCols(), out.L.NumCols());
	for (int i{ 0 }; i < LExpected.NumRows(); i++) {
		for (int j{ 0 }; j < LExpected.NumCols(); j++) {
			EXPECT_LE(relError(LExpected(i, j), out.L(i, j)), rTol) <<
				"In L(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, three_wires_ribbon_generalized_capacitance)
{
	// Three wires ribbon open problem. 
	// Comparison with Clayton Paul's book:  
	// Analysis of multiconductor transmision lines. 2007.
	// Sec. 5.2.3, p. 187.

	const std::string CASE{ "three_wires_ribbon" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };

	auto dr = Driver::loadFromAdaptedFile(fn);

	double gCExpectedData[9] = {
	  26.2148,  -18.0249,  -5.03325,
	 -18.0249,   37.8189, -18.0249,
	  -5.03325, -18.0249,  26.2148
	};
	mfem::DenseMatrix gCExpected(3, 3);
	gCExpected.UseExternalData(gCExpectedData, 3, 3);
	gCExpected *= 1e-12;

	auto gC = dr.getGeneralizedCMatrix();
	gC *= EPSILON0_SI;

	const double rTol{ 0.015 };
	ASSERT_EQ(gCExpected.NumRows(), gC.NumRows());
	ASSERT_EQ(gCExpected.NumCols(), gC.NumCols());
	for (int i{ 0 }; i < gCExpected.NumRows(); i++) {
		for (int j{ 0 }; j < gCExpected.NumCols(); j++) {
			EXPECT_LE(relError(gCExpected(i, j), gC(i, j)), rTol) <<
				"In gC(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, three_wires_ribbon_floating_potentials)
{
	// Three wires ribbon open problem. 
	// Comparison with Clayton Paul's book:  
	// Analysis of multiconductor transmision lines. 2007.
	// Sec. 5.2.3, p. 187.

	const std::string CASE{ "three_wires_ribbon" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };

	auto dr = Driver::loadFromAdaptedFile(fn);
	auto fp = dr.getFloatingPotentials(1, false);

	// Solves problem and checks that charge is zero in the floating conductor.
	auto sP = dr.getElectrostaticSolvedProblem();
	dr.loadFloatingPotentials(sP, fp);
	auto& s = sP->solver;
	// For debugging.
	// ParaViewDataCollection pd(outFolder() + CASE + "_floating", s.getMesh());
	// s.writeParaViewFields(pd);
	auto materials = &dr.getModel().getMaterials();
	auto Q0 = s->getChargeInBoundary(materials->getConductorWithId(0)->getAttribute());
	auto Q1 = s->getChargeInBoundary(materials->getConductorWithId(1)->getAttribute());
	auto Q2 = s->getChargeInBoundary(materials->getConductorWithId(2)->getAttribute());

	auto openBdr = materials->getOpenBoundaries().front();
	auto Qb = s->getChargeInBoundary(openBdr->getAttribute());

	const double aTol{ 1e-3 };
	EXPECT_NEAR(0.0, Q0, aTol);
	EXPECT_NEAR(0.0, Q2, aTol);
	EXPECT_NEAR(0.0, Q0 + Q1 + Q2 + Qb, aTol);

}

TEST_F(DriverTest, nested_coax)
{
	auto out{ Driver::loadFromAdaptedFile(inputCase("nested_coax")).getPULMTL() };


	auto C01{ EPSILON0_SI * 2.0 * M_PI / log(8.0 / 5.6) };
	auto C12{ EPSILON0_SI * 2.0 * M_PI / log(4.8 / 2.0) };

	double CExpectedData[4] = {
		  C01 + C12, -C12,
		 -C12,      C12
	};
	mfem::DenseMatrix CExpected(2, 2);
	CExpected.UseExternalData(CExpectedData, 2, 2);

	const double rTol{ 0.10 };

	ASSERT_EQ(CExpected.NumRows(), out.C.NumRows());
	ASSERT_EQ(CExpected.NumCols(), out.C.NumCols());
	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
			EXPECT_LE(relError(CExpected(i, j), out.C(i, j)), rTol) <<
				"In C(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, DISABLED_nested_coax_by_domains) // Fails when enforcing conductors to start in zero and be consecutive.
{
	auto out{ Driver::loadFromAdaptedFile(inputCase("nested_coax")).getPULMTLByDomains() };

	auto C01{ EPSILON0_SI * 2.0 * M_PI / log(8.0 / 5.6) };
	auto C12{ EPSILON0_SI * 2.0 * M_PI / log(4.8 / 2.0) };

	const double rTol{ 0.10 };

	EXPECT_EQ(2, out.domainTree.verticesSize());
	std::vector<std::pair<int, int>> domainConnections{ std::make_pair(0,1) };
	EXPECT_EQ(domainConnections, out.domainTree.getEdgesAsPairs());

	ASSERT_EQ(2, out.domainToPUL.size());

	ASSERT_EQ(1, out.domainToPUL.at(0).C.NumRows());
	EXPECT_LE(relError(C01, out.domainToPUL.at(0).C(0, 0)), rTol);

	ASSERT_EQ(1, out.domainToPUL.at(1).C.NumRows());
	EXPECT_LE(relError(C12, out.domainToPUL.at(1).C(0, 0)), rTol);
}

TEST_F(DriverTest, agrawal1981)
{
	// Agrawal, Ashok K. and Price, Harold J.
	// Experimental Characterization of Partially Degenerate Three-Conductor
	// Transmission Lines in the Time Domain. IEEE-TEMC. 1981.

	const std::string CASE{ "agrawal1981" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };

	// P.U.L. Capacitances obtained from eigenvectors.
	const int nConductors = 3;
	double CExpectedData[nConductors * nConductors] = {
		  74.54, -34.57, -34.2,
		 -34.63,  73.87, -33.96,
		 -34.29, -34.0,   73.41
	};
	mfem::DenseMatrix CExpected(nConductors, nConductors);
	CExpected.UseExternalData(CExpectedData, nConductors, nConductors);
	CExpected *= 1e-12;

	auto out{ Driver::loadFromAdaptedFile(fn).getPULMTL() };

	const double rTol{ 0.10 };

	ASSERT_EQ(CExpected.NumRows(), out.C.NumRows());
	ASSERT_EQ(CExpected.NumCols(), out.C.NumCols());
	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
			EXPECT_LE(relError(CExpected(i, j), out.C(i, j)), rTol) <<
				"In C(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, lansink2024_floating_potentials)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.

	const std::string CASE{ "lansink2024" };

	auto dr{ Driver::loadFromAdaptedFile(casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json") };

	auto fp = dr.getFloatingPotentials(0, false);
	auto inCell{ dr.getInCellPotentials() };

	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	SolverInputs p;
	p.dirichletBoundaries = {
		{
			{1, fp.at(0)}, // Conductor 0 floating potential.
			{2, fp.at(1)}, // Conductor 1 prescribed potential.
		}
	};
	p.openBoundaries = { 3 };

	tulip::Solver s{ m, p };
	s.Solve();

	auto Q0 = s.getChargeInBoundary(1);
	auto Q1 = s.getChargeInBoundary(2);
	auto Qb = s.getChargeInBoundary(3);

	// For debugging.
	ParaViewDataCollection pd(outFolder() + CASE + "_floating", s.getMesh());
	s.writeParaViewFields(pd);

	// Expectations.
	const double aTol{ 1e-3 };
	EXPECT_NEAR(0.0, Q1, aTol);
	EXPECT_NEAR(0.0, Q0 + Q1 + Qb, aTol);

	const double a0 = inCell.electric.at(0).ab[0].first;
	EXPECT_NEAR(Q0, a0, 1e-4);
}

TEST_F(DriverTest, lansink2024_fdtd_in_cell_parameters_around_conductor_1)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.

	const std::string CASE{ "lansink2024_fdtd_cell" };

	auto inCell{
		Driver::loadFromAdaptedFile(
			casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json"
		).getInCellPotentials()
	};

	const double rTol = 0.03;

	// In this test case inner region coincides with fdtd-cell.
	// In-cell capacitances.
	{
		auto computedC00 = inCell.getCapacitanceUsingInnerRegion(0, 0);
		auto expectedC00 = 14.08e-12; // C11 for floating in paper. Table 1.
		EXPECT_NEAR(0.0, relError(expectedC00, computedC00), rTol);
	}

	{
		auto computedC01 = inCell.getCapacitanceUsingInnerRegion(0, 1);
		auto expectedC01 = 43.99e-12; // C12 for floating in paper. Table 1.
		EXPECT_NEAR(0.0, relError(expectedC01, computedC01), rTol);
	}

	// In-cell inductances
	{
		auto computedL00 = inCell.getInductanceUsingInnerRegion(0, 0);
		auto expectedL00 = 791e-9; // L11 for floating in paper. Table 1.
		EXPECT_NEAR(0.0, relError(expectedL00, computedL00), rTol);
	}

	{
		auto computedL01 = inCell.getInductanceUsingInnerRegion(0, 1);
		auto expectedL01 = 253e-9; // L12 for floating in paper. Table 1.
		EXPECT_NEAR(0.0, relError(expectedL01, computedL01), rTol);
	}
}

TEST_F(DriverTest, lansink2024_two_wires_using_multipolar_expansion)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.

	const std::string CASE{ "lansink2024" };

	auto inCell{
		Driver::loadFromAdaptedFile(
			casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json"
		).getInCellPotentials()
	};


	const double rTol = 0.06;

	Box fdtdCell0{ {-0.110, -0.100}, {0.090, 0.100} };

	auto computedC00 = inCell.getCapacitanceOnBox(0, 0, fdtdCell0);
	auto expectedC00 = 14.08e-12; // C11 for floating in paper. Table 1.
	EXPECT_NEAR(0.0, relError(expectedC00, computedC00), rTol);

	auto computedC01 = inCell.getCapacitanceOnBox(0, 1, fdtdCell0);
	auto expectedC01 = 43.99e-12; // C12 for floating in paper. Table 1.
	EXPECT_NEAR(0.0, relError(expectedC01, computedC01), rTol);

	auto computedL00 = inCell.getInductanceOnBox(0, 0, fdtdCell0);
	auto expectedL00 = 791e-9; // L11 for floating in paper. Table 1.
	EXPECT_NEAR(0.0, relError(expectedL00, computedL00), rTol);

	auto computedL01 = inCell.getInductanceOnBox(0, 1, fdtdCell0);
	auto expectedL01 = 253e-9; // L12 for floating in paper. Table 1.
	EXPECT_NEAR(0.0, relError(expectedL01, computedL01), rTol);

	Box fdtdCell1{ {-0.090, -0.100}, {0.110, 0.100} };

	auto computedC10 = inCell.getCapacitanceOnBox(1, 0, fdtdCell1);
	auto expectedC10 = 44.31e-12; // C21 for floating in paper. Table 1.
	EXPECT_NEAR(0.0, relError(expectedC10, computedC10), rTol);

	auto computedC11 = inCell.getCapacitanceOnBox(1, 1, fdtdCell1);
	auto expectedC11 = 28.79e-12; // C22 for floating in paper. Table 1.
	EXPECT_NEAR(0.0, relError(expectedC11, computedC11), rTol);

	auto computedL10 = inCell.getInductanceOnBox(1, 0, fdtdCell1);
	auto expectedL10 = 251e-9; // L21 for floating in paper. Table 1.
	EXPECT_NEAR(0.0, relError(expectedL10, computedL10), rTol);

	auto computedL11 = inCell.getInductanceOnBox(1, 1, fdtdCell1);
	auto expectedL11 = 387e-9; // L22 for floating in paper. Table 1.
	EXPECT_NEAR(0.0, relError(expectedL11, computedL11), rTol);
}

TEST_F(DriverTest, lansink2024_fdtd_cell_shifted_and_centered)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.

	const std::string CASE{ "lansink2024" };

	auto inCell{
		Driver::loadFromAdaptedFile(
			casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json"
		).getInCellPotentials()
	};

	Box fdtdCellCentered{ {-0.100, -0.100}, {0.100, 0.100} };

	{
		Box fdtdCellShifted{ {-0.110, -0.100}, {0.090, 0.100} };
		auto computedC_shifted = inCell.getCapacitanceOnBox(0, 0, fdtdCellShifted);
		auto computedC_centered = inCell.getCapacitanceOnBox(0, 0, fdtdCellCentered);
		auto err = relError(computedC_shifted, computedC_centered);
		EXPECT_TRUE(err < 1e-4);
	}
	{
		Box fdtdCellShifted{ {-0.110, -0.100}, {0.090, 0.100} };
		auto computedC_shifted = inCell.getCapacitanceOnBox(0, 1, fdtdCellShifted);
		auto computedC_centered = inCell.getCapacitanceOnBox(0, 1, fdtdCellCentered);
		auto err = relError(computedC_shifted, computedC_centered);
		EXPECT_TRUE(err < 1e-2);
	}
	{
		Box fdtdCellShifted{ { -0.090, -0.100 }, { 0.110, 0.100 } };
		auto computedC_shifted = inCell.getCapacitanceOnBox(1, 0, fdtdCellShifted);
		auto computedC_centered = inCell.getCapacitanceOnBox(1, 0, fdtdCellCentered);
		auto err = relError(computedC_shifted, computedC_centered);
		EXPECT_TRUE(err < 1e-2);
	}
	{
		Box fdtdCellShifted{ { -0.090, -0.100 }, { 0.110, 0.100 } };
		auto computedC_shifted = inCell.getCapacitanceOnBox(1, 1, fdtdCellShifted);
		auto computedC_centered = inCell.getCapacitanceOnBox(1, 1, fdtdCellCentered);
		auto err = relError(computedC_shifted, computedC_centered);
		EXPECT_TRUE(err < 1e-2);
	}
}

TEST_F(DriverTest, lansink2024_single_wire_in_cell_parameters)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.
	// VALUES IN TABLE 3 ARE WRONG AND HAVE BEEN CORRECTED.

	const std::string CASE{ "lansink2024_single_wire" };

	auto inCell{
		Driver::loadFromAdaptedFile(
			casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json"
		).getInCellPotentials()
	};

	const double rTol = 0.05;

	// In this test case inner region coincides with fdtd-cell.
	// In-cell capacitances.
	{
		auto computedC00 = inCell.getCapacitanceUsingInnerRegion(0, 0);
		auto expectedC00 = 49.11e-12; // C11 with insulation. Table 3. 
		// Paper has a mistake, this is the correct value.
		EXPECT_NEAR(0.0, relError(expectedC00, computedC00), rTol);
	}

	// In-cell inductances
	{
		auto computedL00 = inCell.getInductanceUsingInnerRegion(0, 0);
		auto expectedL00 = 320e-9; // L11 with insulation. Table 3. 
		// Paper has a mistake, this is the correct value.
		EXPECT_NEAR(0.0, relError(expectedL00, computedL00), rTol);
	}

}

TEST_F(DriverTest, lansink2024_single_wire_multipolar_in_cell_parameters)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.
	// VALUES IN TABLE 3 ARE WRONG AND HAVE BEEN CORRECTED.

	const std::string CASE{ "lansink2024_single_wire_multipolar" };

	auto inCell{
		Driver::loadFromAdaptedFile(
			casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json"
		).getInCellPotentials()
	};

	Box fdtdCell{ {-0.0075, -0.0075}, {0.0075, 0.0075} };

	const double rTol = 0.06;

	// In this test case inner region coincides with fdtd-cell.
	// In-cell capacitances.
	{
		auto computedC00 = inCell.getCapacitanceOnBox(0, 0, fdtdCell);
		auto expectedC00 = 49.11e-12; // C11 with insulation. Table 3. 
		// Paper has a mistake, this is the correct value.
		EXPECT_NEAR(0.0, relError(expectedC00, computedC00), rTol);
	}

	// In-cell inductances
	{
		auto computedL00 = inCell.getInductanceOnBox(0, 0, fdtdCell);
		auto expectedL00 = 320e-9; // L11 with insulation. Table 3. 
		// Paper has a mistake, this is the correct value.
		EXPECT_NEAR(0.0, relError(expectedL00, computedL00), rTol);
	}

	// Check that multipolar expansion for the bare wire produces a 1 V at the boundary.
	auto a0 = inCell.magnetic.at(0).ab[0].first;
	auto Va = a0 / (2 * M_PI) * log(1.0 / 1e-3);
	EXPECT_NEAR(1.0, Va, 1e-3);

	saveToJSONFile(inCell.toJSON(),
		"lansink2024_single_wire_multipolar.inCellPotentials.out.json");
}

TEST_F(DriverTest, getCFromGeneralizedC_two_wires_open)
{
	const std::string CASE{ "two_wires_open" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };

	auto driver = Driver::loadFromAdaptedFile(fn);
	auto gC{ driver.getGeneralizedCMatrix() };

	double d = 50;
	double rw1 = 2;
	double rw2 = 2;

	double CExpected{
		2 * M_PI * EPSILON0_NATURAL /
		std::acosh((d * d - rw1 * rw1 - rw2 * rw2) / (2 * rw1 * rw2))
	};

	auto C = driver.getCFromGeneralizedC(gC, Model::Openness::open);

	const double rTol{ 0.0025 };
	ASSERT_EQ(1, C.NumRows());
	ASSERT_EQ(1, C.NumCols());
	EXPECT_LE(relError(CExpected, C(0, 0)), rTol) <<
		"In C(" << 0 << ", " << 0 << ")";
}

TEST_F(DriverTest, getCFromGeneralizedC_three_wires)
{
	// Three wires ribbon open problem. 
	// Comparison with Clayton Paul's book:  
	// Analysis of multiconductor transmision lines. 2007.
	// Sec. 5.2.3, p. 187.
	const std::string CASE{ "two_wires_open" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };
	auto driver = Driver::loadFromAdaptedFile(fn);

	double gCData[9] = {
	  26.2148,  -18.0249,  -5.03325,
	 -18.0249,   37.8189, -18.0249,
	  -5.03325, -18.0249,  26.2148
	};
	mfem::DenseMatrix gC(3, 3);
	gC.UseExternalData(gCData, 3, 3);
	gC *= 1e-12;


	double CExpectedData[4] = {
		  37.432, -18.716,
		 -18.716,  24.982
	};
	mfem::DenseMatrix CExpected(2, 2);
	CExpected.UseExternalData(CExpectedData, 2, 2);
	CExpected *= 1e-12;




	auto C = driver.getCFromGeneralizedC(gC, Model::Openness::open);

	const double rTol{ 0.001 };
	ASSERT_EQ(CExpected.NumRows(), C.NumRows());
	ASSERT_EQ(CExpected.NumCols(), C.NumCols());
	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
			EXPECT_LE(relError(CExpected(i, j), C(i, j)), rTol) <<
				"In C(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, lansink2024_small_one_centered_fdtd_cell_vs_multipolar)
{
	// In-cell capacitances centered in conductor 0
	// Using meshed FDTD cell.
	InCellPotentials fdtdCellPotentials;
	{
		const std::string CASE{ "lansink2024_small_one_centered_fdtd_cell" };
		fdtdCellPotentials = Driver::loadFromAdaptedFile(
			casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json").getInCellPotentials();
	}
	auto fdtdCellComputedC00 = fdtdCellPotentials.getCapacitanceUsingInnerRegion(0, 0);
	auto fdtdCellComputedC01 = fdtdCellPotentials.getCapacitanceUsingInnerRegion(0, 1);

	// Using multipolar expansion.
	InCellPotentials multipolarPotentials;
	{
		const std::string CASE{ "lansink2024_small_one_centered" };
		multipolarPotentials = Driver::loadFromAdaptedFile(
			casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json").getInCellPotentials();
	}
	Box fdtdCell{ {-0.1, -0.1}, {0.1, 0.1} };
	auto multipolarComputedC00 = multipolarPotentials.getCapacitanceOnBox(0, 0, fdtdCell);
	auto multipolarComputedC01 = multipolarPotentials.getCapacitanceOnBox(0, 1, fdtdCell);

	// Compares results
	EXPECT_NEAR(fdtdCellComputedC00, multipolarComputedC00, 0.2e-12);
	EXPECT_NEAR(fdtdCellComputedC01, multipolarComputedC01, 0.2e-12);
	EXPECT_LE(relError(fdtdCellComputedC00, multipolarComputedC00), 0.005);
	EXPECT_LE(relError(fdtdCellComputedC01, multipolarComputedC01), 0.005);

	saveToJSONFile(multipolarPotentials.toJSON(),
		"lansink2024_small_one_centered.inCellPotentials.out.json");
}

TEST_F(DriverTest, lansink2024_large_one_centered_fdtd_cell)
{
	// In-cell capacitances centered in conductor 0
	// Using meshed FDTD cell.
	InCellPotentials fdtdCellPotentials;
	{
		const std::string CASE{ "lansink2024_large_one_centered_fdtd_cell" };
		fdtdCellPotentials = Driver::loadFromAdaptedFile(
			casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json").getInCellPotentials();
	}
	auto fdtdCellComputedC10 = fdtdCellPotentials.getCapacitanceUsingInnerRegion(1, 0);

	// Compares results
	EXPECT_NEAR(fdtdCellComputedC10, 43.8646e-12, 1e-12); // Computed with BEM.
}

TEST_F(DriverTest, lansink2024_small_one_centered_different_integration_centers)
{
	// In-cell capacitances centered in conductor 0
	InCellPotentials multipolarPotentials;
	const std::string CASE{ "lansink2024_small_one_centered" };

	multipolarPotentials = Driver::loadFromAdaptedFile(
		casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json").getInCellPotentials();

	mfem::DenseMatrix geometricC(2, 2);
	Box fdtdCellCenteredOnConductor0{ {-0.1, -0.1}, {0.1, 0.1} };
	{
		geometricC(0, 0) = multipolarPotentials.getCapacitanceOnBox(0, 0, fdtdCellCenteredOnConductor0);
		geometricC(0, 1) = multipolarPotentials.getCapacitanceOnBox(0, 1, fdtdCellCenteredOnConductor0);
	}
	{
		Box fdtdCellCenteredOnConductor1{ {-0.12, -0.1}, {0.08, 0.1} };
		geometricC(1, 0) = multipolarPotentials.getCapacitanceOnBox(1, 0, fdtdCellCenteredOnConductor1);
		geometricC(1, 1) = multipolarPotentials.getCapacitanceOnBox(1, 1, fdtdCellCenteredOnConductor1);
	}

	mfem::DenseMatrix chargeCenteredC(2, 2);
	for (int i = 0; i < 2; i++) {
		Box chargeCenteredCell{ fdtdCellCenteredOnConductor0 };
		chargeCenteredCell.displace(multipolarPotentials.electric.at(i).expansionCenter);
		for (int j = 0; j < 2; j++) {
			chargeCenteredC(i, j) = multipolarPotentials.getCapacitanceOnBox(i, j, chargeCenteredCell);
		}
	}

	// Compares results
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			EXPECT_NEAR(geometricC(i, j), chargeCenteredC(i, j), 0.1e-11);
		}
	}

}

TEST_F(DriverTest, DISABLED_realistic_case_with_dielectrics_fdtd_cell)
{
	const std::string CASE{ "realistic_case_with_dielectrics_fdtd_cell" };
	auto dr = Driver::loadFromAdaptedFile(
		casesFolder() + "realistic_case_with_dielectrics/" + CASE + ".tulip.adapted.json");
	auto fdtdCellPotentials = dr.getInCellPotentials();
	auto fdtdCellComputed_C_0 = fdtdCellPotentials.getCapacitanceUsingInnerRegion(0, 0);
	auto fdtdCellComputed_C_16 = fdtdCellPotentials.getCapacitanceUsingInnerRegion(0, 16);
	auto fdtdCellComputed_C_25 = fdtdCellPotentials.getCapacitanceUsingInnerRegion(0, 25);
	auto fdtdCellComputed_C_30 = fdtdCellPotentials.getCapacitanceUsingInnerRegion(0, 30);

	auto fdtdCellComputed_L_0 = fdtdCellPotentials.getInductanceUsingInnerRegion(0, 0);
	auto fdtdCellComputed_L_16 = fdtdCellPotentials.getInductanceUsingInnerRegion(0, 16);
	auto fdtdCellComputed_L_25 = fdtdCellPotentials.getInductanceUsingInnerRegion(0, 25);
	auto fdtdCellComputed_L_30 = fdtdCellPotentials.getInductanceUsingInnerRegion(0, 30);

	auto rTol = 0.001;
	EXPECT_LE(relError(fdtdCellComputed_C_0, 4.0911726228481947e-11), rTol);
	EXPECT_LE(relError(fdtdCellComputed_C_16, 1.2547925523607968e-10), rTol);
	EXPECT_LE(relError(fdtdCellComputed_C_25, 5.9060595987059621e-11), rTol);
	EXPECT_LE(relError(fdtdCellComputed_C_30, 1.7797154919313720e-10), rTol);
	EXPECT_LE(relError(fdtdCellComputed_L_0, 2.9989786293920517e-07), rTol);
	EXPECT_LE(relError(fdtdCellComputed_L_16, 8.7458344767957649e-08), rTol);
	EXPECT_LE(relError(fdtdCellComputed_L_25, 2.0233814651758583e-07), rTol);
	EXPECT_LE(relError(fdtdCellComputed_L_30, 5.9647510363871327e-08), rTol);
}

TEST_F(DriverTest, DISABLED_realistic_case_with_dielectrics_multipolar)
{
	
	const std::string CASE{ "realistic_case_with_dielectrics_inner_region" };
	auto dr = Driver::loadFromAdaptedFile(
		casesFolder() + "realistic_case_with_dielectrics/" + CASE + ".tulip.adapted.json");
	auto mP = dr.getInCellPotentials();
	
	Box fdtdCellCenteredOnConductor0{ {-0.016209, -0.009066}, {0.013791, 0.020934} };
	auto mPComputedC_0 = mP.getCapacitanceOnBox(0, 0, fdtdCellCenteredOnConductor0);
	auto mPComputedC_16 = mP.getCapacitanceOnBox(0, 16, fdtdCellCenteredOnConductor0);
	auto mPComputedC_25 = mP.getCapacitanceOnBox(0, 25, fdtdCellCenteredOnConductor0);
	auto mPComputedC_30 = mP.getCapacitanceOnBox(0, 30, fdtdCellCenteredOnConductor0);

	auto mPComputedL_0 = mP.getInductanceOnBox(0, 0, fdtdCellCenteredOnConductor0);
	auto mPComputedL_16 = mP.getInductanceOnBox(0, 16, fdtdCellCenteredOnConductor0);
	auto mPComputedL_25 = mP.getInductanceOnBox(0, 25, fdtdCellCenteredOnConductor0);
	auto mPComputedL_30 = mP.getInductanceOnBox(0, 30, fdtdCellCenteredOnConductor0);

	// Compare with C and L computed using the fdtd cell.
	const double rTol = 0.015;
	const double fdtdCellComputed_C_0 = 4.0911726228481947e-11;
	const double fdtdCellComputed_C_16 = 1.2547925523607968e-10;
	const double fdtdCellComputed_C_25 = 5.9060595987059621e-11;
	const double fdtdCellComputed_C_30 = 1.7797154919313720e-10;
	const double fdtdCellComputed_L_0 = 2.9989786293920517e-07;
	const double fdtdCellComputed_L_16 = 8.7458344767957649e-08;
	const double fdtdCellComputed_L_25 = 2.0233814651758583e-07;
	const double fdtdCellComputed_L_30 = 5.9647510363871327e-08;

	EXPECT_LE(relError(fdtdCellComputed_C_0, mPComputedC_0), rTol) << "C(0,0) mismatch";
	EXPECT_LE(relError(fdtdCellComputed_C_16, mPComputedC_16), rTol) << "C(0,16) mismatch";
	EXPECT_LE(relError(fdtdCellComputed_C_25, mPComputedC_25), rTol) << "C(0,25) mismatch";
	EXPECT_LE(relError(fdtdCellComputed_C_30, mPComputedC_30), rTol) << "C(0,30) mismatch";
	EXPECT_LE(relError(fdtdCellComputed_L_0, mPComputedL_0), rTol) << "L(0,0) mismatch";
	EXPECT_LE(relError(fdtdCellComputed_L_16, mPComputedL_16), rTol) << "L(0,16) mismatch";
	EXPECT_LE(relError(fdtdCellComputed_L_25, mPComputedL_25), rTol) << "L(0,25) mismatch";
	EXPECT_LE(relError(fdtdCellComputed_L_30, mPComputedL_30), rTol) << "L(0,30) mismatch";


	// // For debugging
	// saveToJSONFile(mP.toJSON(),
	// 	"realistic_case_with_dielectrics.inCellPotentials.out.json");
}
