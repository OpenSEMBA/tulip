#include <gtest/gtest.h>

#include "TestUtils.h"

#include "Results.h"

using namespace pulmtln;

using json = nlohmann::json;

class ResultsTest : public ::testing::Test {
};

TEST_F(ResultsTest, PULParameters_serialization_to_JSON)
{
	PULParameters p;
	p.C = mfem::DenseMatrix(2, 2);
	p.C(0, 0) = 1.0; p.C(0, 1) = 2.0;
	p.C(1, 0) = 3.0; p.C(1, 1) = 4.0;

	p.L = mfem::DenseMatrix(2, 2);
	p.L(0, 0) = 1.0; p.L(0, 1) = 2.0;
	p.L(1, 0) = 3.0; p.L(1, 1) = 4.0;

	json j = p.toJSON();
	
	PULParameters r(j);

	EXPECT_EQ(p, r);
}

TEST_F(ResultsTest, PULParameters_serialization_to_JSON_with_generalized)
{
	PULParameters p;
	p.C = mfem::DenseMatrix(2, 2);
	p.C(0, 0) = 1.0; p.C(0, 1) = 2.0;
	p.C(1, 0) = 3.0; p.C(1, 1) = 4.0;

	p.L = mfem::DenseMatrix(2, 2);
	p.L(0, 0) = 1.0; p.L(0, 1) = 2.0;
	p.L(1, 0) = 3.0; p.L(1, 1) = 4.0;

	p.gC = mfem::DenseMatrix(3, 3);
	p.gC(0, 0) = 9.0; p.gC(0, 1) = 8.0; p.gC(0, 2) = 7.0;
	p.gC(1, 0) = 6.0; p.gC(1, 1) = 5.0; p.gC(1, 2) = 4.0;
	p.gC(2, 0) = 3.0; p.gC(2, 1) = 2.0; p.gC(2, 2) = 1.0;

	p.gL = mfem::DenseMatrix(3, 3);
	p.gL(0, 0) = 1.0; p.gL(0, 1) = 2.0; p.gL(0, 2) = 3.0;
	p.gL(1, 0) = 4.0; p.gL(1, 1) = 5.0; p.gL(1, 2) = 6.0;
	p.gL(2, 0) = 7.0; p.gL(2, 1) = 8.0; p.gL(2, 2) = 9.0;

	json j = p.toJSON();

	EXPECT_TRUE(j.contains("gC"));
	EXPECT_TRUE(j.contains("gL"));

	PULParameters r(j);

	EXPECT_EQ(p, r);
}

TEST_F(ResultsTest, PULParameters_serialization_backward_compatibility_no_generalized)
{
	// JSON without gC/gL should still deserialize correctly (backward compatibility).
	PULParameters p;
	p.C = mfem::DenseMatrix(1, 1);
	p.C(0, 0) = 5.0;
	p.L = mfem::DenseMatrix(1, 1);
	p.L(0, 0) = 3.0;

	json j = p.toJSON();

	EXPECT_FALSE(j.contains("gC"));
	EXPECT_FALSE(j.contains("gL"));

	PULParameters r(j);

	EXPECT_EQ(0, r.gC.NumRows());
	EXPECT_EQ(0, r.gL.NumRows());
	EXPECT_EQ(p.C(0, 0), r.C(0, 0));
	EXPECT_EQ(p.L(0, 0), r.L(0, 0));
}

TEST_F(ResultsTest, InCellPotentials_serialization_to_JSON)
{
	// Dummy values
	InCellPotentials p;
	p.innerRegionBox = Box{ { 0.0, 0.0 }, { 1.0, 1.0 } };
	
	p.electric[0].innerRegionAveragePotential = 1.0;
	p.electric[0].expansionCenter = { 0.5, 0.5 };
	p.electric[0].ab = {
		{1.0, 0.0}, {0.5, 0.5}, {0.25, 0.25}
	};
	p.electric[0].conductorPotentials[0] = 1.0;
	p.electric[0].conductorPotentials[1] = 2.0;

	p.magnetic[0].innerRegionAveragePotential = 2.0;
	p.magnetic[0].expansionCenter = { 0.5, 0.5 };
	p.magnetic[0].ab = {
		{2.0, 0.0}, {1.0, 1.0}, {0.5, 0.5}
	};
	p.magnetic[0].conductorPotentials[0] = 1.0;
	p.magnetic[0].conductorPotentials[1] = 2.0;

	json j = p.toJSON();
	
	InCellPotentials r(j);

	EXPECT_EQ(p, r);

}

TEST_F(ResultsTest, InCellPotentials_serialization_to_JSON_2)
{
	// Dummy values
	InCellPotentials p;
	p.innerRegionBox = Box{ { 0.0, 0.0 }, { 1.0, 1.0 } };
	
	// Electric
	p.electric[0].innerRegionAveragePotential = 1.0;
	p.electric[0].expansionCenter = { 0.5, 0.5 };
	p.electric[0].ab = { {1.0, 0.0}, {0.5, 0.5}, {0.25, 0.25} };
	p.electric[0].conductorPotentials[0] = 1.0;
	p.electric[0].conductorPotentials[1] = 2.0;

	p.electric[1].innerRegionAveragePotential = 1.0;
	p.electric[1].expansionCenter = { 0.5, 0.5 };
	p.electric[1].ab = {{1.0, 0.0}, {0.5, 0.5}, {0.25, 0.25} };
	p.electric[1].conductorPotentials[0] = 1.0;
	p.electric[1].conductorPotentials[1] = 2.0;

	// Magnetic
	p.magnetic[0].innerRegionAveragePotential = 2.0;
	p.magnetic[0].expansionCenter = { 0.5, 0.5 };
	p.magnetic[0].ab = { {2.0, 0.0}, {1.0, 1.0}, {0.5, 0.5} };
	p.magnetic[0].conductorPotentials[0] = 1.0;
	p.magnetic[0].conductorPotentials[1] = 2.0;

	p.magnetic[1].innerRegionAveragePotential = 2.0;
	p.magnetic[1].expansionCenter = { 0.5, 0.5 };
	p.magnetic[1].ab = { {2.0, 0.0}, {1.0, 1.0}, {0.5, 0.5} };
	p.magnetic[1].conductorPotentials[0] = 1.0;
	p.magnetic[1].conductorPotentials[1] = 2.0;

	json j = p.toJSON();
	InCellPotentials r(j);
	EXPECT_EQ(p, r);

}