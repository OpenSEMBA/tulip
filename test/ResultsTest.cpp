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