#include <gtest/gtest.h>

#include "Model.h"
#include "AdaptedInputParser.h"
#include "TestUtils.h"

using namespace tulip;

class ModelTest : public ::testing::Test {};

TEST_F(ModelTest, empty_coax_is_closed)
{
	EXPECT_EQ(
		Parser{ inputCase("empty_coax") }.readModel().determineOpenness(),
		Model::Openness::closed
	);
}

TEST_F(ModelTest, three_wires_ribbon_is_open)
{
	EXPECT_EQ(
		Parser{ inputCase("three_wires_ribbon") }.readModel().determineOpenness(),
		Model::Openness::open
	);
}

TEST_F(ModelTest, agrawal1981_conductors_in_mesh)
{
	auto m{ Parser{ inputCase("agrawal1981") }.readModel() };
	
	ASSERT_EQ(4, m.getMaterials().getConductors().size());
}

TEST_F(ModelTest, lansink2024_fdtd_cell_material_areas)
{
	const std::string CASE{ "lansink2024_fdtd_cell" };
	auto model{ Parser{ inputCase(CASE) }.readModel() };
	const double rTol{ 1e-8 };

	// Conductor areas (enclosed area computed from boundary integrals)
	Conductor conductor0{ 1, 0 };
	Conductor conductor1{ 2, 1 };
	EXPECT_NEAR(0.0, relError(3.14159253264265e-06, model.getAreaOfMaterial(&conductor0)), rTol);
	EXPECT_NEAR(0.0, relError(3.14159253264264e-04, model.getAreaOfMaterial(&conductor1)), rTol);

	// Dielectric domain areas
	Dielectric dielectric4{ 4 };
	Dielectric dielectric5{ 5 };
	EXPECT_NEAR(0.0, relError(0.0617876225248815, model.getAreaOfMaterial(&dielectric4)), rTol);
	EXPECT_NEAR(0.0, relError(0.0396826990779316, model.getAreaOfMaterial(&dielectric5)), rTol);
}

