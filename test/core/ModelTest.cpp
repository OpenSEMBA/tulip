#include <gtest/gtest.h>

#include "Model.h"
#include "Parser.h"
#include "TestUtils.h"

using namespace pulmtln;

class ModelTest : public ::testing::Test {};

TEST_F(ModelTest, empty_coax_is_closed)
{
	EXPECT_EQ(
		Parser{ inputCase("empty_coax") }.readModel().determineOpenness(),
		Model::Openness::closed
	);
}

TEST_F(ModelTest, agrawal1981_is_semiopen)
{
	EXPECT_EQ(
		Parser{ inputCase("agrawal1981") }.readModel().determineOpenness(),
		Model::Openness::semiopen
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
	
	ASSERT_EQ(4, m.getMaterials().pecs.size());
}

TEST_F(ModelTest, lansink2024_fdtd_cell_areas)
{
	auto m{ Parser{ inputCase("lansink2024_fdtd_cell") }.readModel() };

	auto c0 = m.getAreaOfMaterial("Conductor_0");
	auto c1 = m.getAreaOfMaterial("Conductor_1");
	auto v0 = m.getAreaOfMaterial("Vacuum_0");

	EXPECT_NEAR(0.04, c0 + c1 + v0, 1e-8);
}

TEST_F(ModelTest, lansink2024_fdtd_cell_bounding_box)
{
	auto m{ Parser{ inputCase("lansink2024_fdtd_cell") }.readModel() };
	auto box = m.getBoundingBoxOfMaterial("Vacuum_0");
	EXPECT_NEAR(-0.11, box.min[0], 1e-8);
	EXPECT_NEAR(-0.10, box.min[1], 1e-8);
	EXPECT_NEAR( 0.09, box.max[0], 1e-8);
	EXPECT_NEAR( 0.10, box.max[1], 1e-8);
}
