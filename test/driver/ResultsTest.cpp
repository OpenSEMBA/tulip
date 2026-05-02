#include <gtest/gtest.h>

#include "TestUtils.h"

#include "Results.h"

using namespace tulip;

using json = nlohmann::json;

class ResultsTest : public ::testing::Test {
};

TEST_F(ResultsTest, PULParameters_serialization_deserialization_to_JSON)
{
	PULParameters p;
	p.R = mfem::Vector(2);
	p.R[0] = 5.0; p.R[1] = 0.0;
	p.C = mfem::DenseMatrix(2, 2);
	p.C(0, 0) = 1.0; p.C(0, 1) = 2.0;
	p.C(1, 0) = 3.0; p.C(1, 1) = 4.0;

	p.L = mfem::DenseMatrix(2, 2);
	p.L(0, 0) = 1.0; p.L(0, 1) = 2.0;
	p.L(1, 0) = 3.0; p.L(1, 1) = 4.0;

	json j = p.toJSON(); // Serialization.
	
	PULParameters r(j); // Deserialization.

	EXPECT_EQ(p, r);
}

TEST_F(ResultsTest, PULParameters_deserialization_defaults_missing_R_to_zero)
{
	json j = {
		{"C", {{1.0, 0.0}, {0.0, 2.0}}},
		{"L", {{3.0, 0.0}, {0.0, 4.0}}}
	};

	PULParameters p(j);

	ASSERT_EQ(2, p.R.Size());
	EXPECT_DOUBLE_EQ(0.0, p.R[0]);
	EXPECT_DOUBLE_EQ(0.0, p.R[1]);
}

TEST_F(ResultsTest, InCellPotentials_serialization_to_JSON)
{
	// Dummy values
	InCellPotentials p;
	p.getInnerRegionBox() = Box{ { 0.0, 0.0 }, { 1.0, 1.0 } };
	
	p.getElectric()[0].innerRegionAveragePotential = 1.0;
	p.getElectric()[0].expansionCenter = { 0.5, 0.5 };
	p.getElectric()[0].ab = {
		{1.0, 0.0}, {0.5, 0.5}, {0.25, 0.25}
	};
	p.getElectric()[0].conductorPotentials[0] = 1.0;
	p.getElectric()[0].conductorPotentials[1] = 2.0;

	p.getMagnetic()[0].innerRegionAveragePotential = 2.0;
	p.getMagnetic()[0].expansionCenter = { 0.5, 0.5 };
	p.getMagnetic()[0].ab = {
		{2.0, 0.0}, {1.0, 1.0}, {0.5, 0.5}
	};
	p.getMagnetic()[0].conductorPotentials[0] = 1.0;
	p.getMagnetic()[0].conductorPotentials[1] = 2.0;

	json j = p.toJSON();
	
	InCellPotentials r(j);

	EXPECT_EQ(p, r);

}

TEST_F(ResultsTest, InCellPotentials_serialization_to_JSON_2)
{
	// Dummy values
	InCellPotentials p;
	p.getInnerRegionBox() = Box{ { 0.0, 0.0 }, { 1.0, 1.0 } };
	
	// Electric
	p.getElectric()[0].innerRegionAveragePotential = 1.0;
	p.getElectric()[0].expansionCenter = { 0.5, 0.5 };
	p.getElectric()[0].ab = { {1.0, 0.0}, {0.5, 0.5}, {0.25, 0.25} };
	p.getElectric()[0].conductorPotentials[0] = 1.0;
	p.getElectric()[0].conductorPotentials[1] = 2.0;

	p.getElectric()[1].innerRegionAveragePotential = 1.0;
	p.getElectric()[1].expansionCenter = { 0.5, 0.5 };
	p.getElectric()[1].ab = {{1.0, 0.0}, {0.5, 0.5}, {0.25, 0.25} };
	p.getElectric()[1].conductorPotentials[0] = 1.0;
	p.getElectric()[1].conductorPotentials[1] = 2.0;

	// Magnetic
	p.getMagnetic()[0].innerRegionAveragePotential = 2.0;
	p.getMagnetic()[0].expansionCenter = { 0.5, 0.5 };
	p.getMagnetic()[0].ab = { {2.0, 0.0}, {1.0, 1.0}, {0.5, 0.5} };
	p.getMagnetic()[0].conductorPotentials[0] = 1.0;
	p.getMagnetic()[0].conductorPotentials[1] = 2.0;

	p.getMagnetic()[1].innerRegionAveragePotential = 2.0;
	p.getMagnetic()[1].expansionCenter = { 0.5, 0.5 };
	p.getMagnetic()[1].ab = { {2.0, 0.0}, {1.0, 1.0}, {0.5, 0.5} };
	p.getMagnetic()[1].conductorPotentials[0] = 1.0;
	p.getMagnetic()[1].conductorPotentials[1] = 2.0;

	json j = p.toJSON();
	InCellPotentials r(j);
	EXPECT_EQ(p, r);

}